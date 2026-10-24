#!/usr/bin/env python3
"""
encode_tf_ref.py — Rebuild a TF-binding-site reference from ENCODE.

Steps (each is a subcommand; `all` runs 1-3 in order):
  1. search    query the ENCODE portal, save the raw search JSON
  2. extract   filter file records -> manifest.tsv, optionally download (+md5 check)
  3. stats     per-TF / per-experiment summary tables (+ optional output-type coverage)
  4. (R)       tf_merge_to_grl.R merges peaks per TF into a GRangesList

Typical use:
  python encode_tf_ref.py all --assembly mm10 --outdir ref_mm10 --download -j 8
  Rscript tf_merge_to_grl.R --manifest ref_mm10/manifest.tsv --outprefix ref_mm10/mm10_TF

Each step can also be imported and called from Python (see the notebook).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from urllib.parse import urlencode, urljoin

import requests

BASE = "https://www.encodeproject.org"
HEADERS = {"accept": "application/json"}

DEFAULT_OUTPUT_TYPE = "conservative IDR thresholded peaks"

# Fields requested from the search endpoint. Asking for explicit fields keeps the
# response small and guarantees every record has the same shape, so no per-file
# re-fetch is needed (the old script fell back to one GET per file).
SEARCH_FIELDS = [
    "accession", "status", "href", "md5sum", "file_size",
    "file_format", "file_type", "output_type", "assembly",
    "dataset", "target.label", "biosample_ontology.term_name",
    "biosample_ontology.classification", "biological_replicates",
    "preferred_default", "date_created", "award.rfa",
    "analyses.pipeline_version",
]

MANIFEST_COLS = [
    "target", "experiment_accession", "file_accession", "biosample",
    "biosample_type", "file_format", "file_type", "output_type", "assembly",
    "biological_replicates", "preferred_default", "date_created",
    "pipeline_version", "file_size", "md5sum", "href", "local_path",
    "download_status",
]


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def session() -> requests.Session:
    s = requests.Session()
    s.headers.update(HEADERS)
    adapter = requests.adapters.HTTPAdapter(
        max_retries=requests.adapters.Retry(
            total=5, backoff_factor=2, status_forcelist=(429, 500, 502, 503, 504)
        ),
        pool_maxsize=32,
    )
    s.mount("https://", adapter)
    return s


def build_search_url(assembly: str, output_types: list[str], rfa: str | None,
                     assay_title: str = "TF ChIP-seq", file_formats: list[str] | None = None,
                     organism: str | None = None, fields: list[str] | None = None) -> str:
    params: list[tuple[str, str]] = [
        ("type", "File"),
        ("assay_title", assay_title),
        ("assembly", assembly),
        ("status", "released"),
    ]
    params += [("output_type", ot) for ot in output_types]
    if rfa:
        params.append(("award.rfa", rfa))
    for ff in file_formats or []:
        params.append(("file_format", ff))
    if organism:
        params.append(("replicate.library.biosample.donor.organism.scientific_name", organism))
    for f in fields or SEARCH_FIELDS:
        params.append(("field", f))
    params += [("limit", "all"), ("format", "json")]
    return f"{BASE}/search/?{urlencode(params)}"


def get_nested(d: dict, dotted: str):
    """get_nested(rec, 'target.label'); lists of dicts -> first non-empty value."""
    cur = d
    for key in dotted.split("."):
        if isinstance(cur, list):
            cur = next((x.get(key) for x in cur if isinstance(x, dict) and x.get(key)), None)
        elif isinstance(cur, dict):
            cur = cur.get(key)
        else:
            return None
        if cur is None:
            return None
    return cur


def md5_of(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def write_tsv(rows: list[dict], path: str, cols: list[str]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_tsv(path: str) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ----------------------------------------------------------------------------
# step 1: search
# ----------------------------------------------------------------------------
def step_search(assembly: str = "mm10", outdir: str = ".", rfa: str | None = "ENCODE4",
                output_types: list[str] | None = None, organism: str | None = None) -> str:
    output_types = output_types or [DEFAULT_OUTPUT_TYPE]
    os.makedirs(outdir, exist_ok=True)
    url = build_search_url(assembly, output_types, rfa, organism=organism)
    log(f"GET {url}")
    r = session().get(url, timeout=600)
    r.raise_for_status()
    js = r.json()
    n = len(js.get("@graph", []))

    out_json = os.path.join(outdir, "search_files.json")
    with open(out_json, "w") as f:
        json.dump(js, f)
    # provenance: what was queried and when (ENCODE content changes over time)
    with open(os.path.join(outdir, "search_meta.json"), "w") as f:
        json.dump({
            "url": url, "n_records": n, "assembly": assembly, "rfa": rfa,
            "output_types": output_types,
            "queried_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }, f, indent=2)
    log(f"{n} file records -> {out_json}")
    return out_json


# ----------------------------------------------------------------------------
# step 2: extract (+ download)
# ----------------------------------------------------------------------------
def record_to_row(rec: dict) -> dict:
    dataset = rec.get("dataset") or ""
    exp_acc = dataset.strip("/").split("/")[1] if dataset.startswith("/experiments/") else ""
    reps = rec.get("biological_replicates") or []
    return {
        "target": get_nested(rec, "target.label") or "",
        "experiment_accession": exp_acc,
        "file_accession": rec.get("accession", ""),
        "biosample": get_nested(rec, "biosample_ontology.term_name") or "",
        "biosample_type": get_nested(rec, "biosample_ontology.classification") or "",
        "file_format": rec.get("file_format", ""),
        "file_type": rec.get("file_type", ""),
        "output_type": rec.get("output_type", ""),
        "assembly": rec.get("assembly", ""),
        "biological_replicates": ",".join(map(str, sorted(reps))),
        "preferred_default": bool(rec.get("preferred_default")),
        "date_created": rec.get("date_created", ""),
        "pipeline_version": get_nested(rec, "analyses.pipeline_version") or "",
        "file_size": rec.get("file_size", ""),
        "md5sum": rec.get("md5sum", ""),
        "href": rec.get("href", ""),
    }


def dedup_rows(rows: list[dict], mode: str) -> list[dict]:
    """
    ENCODE sometimes has >1 conservative IDR file per experiment (re-analyses with
    newer pipeline versions, or GRCh38/mm10 lift-overs). Options:
      none       keep everything
      preferred  per (experiment, format): keep preferred_default files if any exist,
                 otherwise keep all of that group
      latest     per (experiment, format): keep the most recently created file
    """
    if mode == "none":
        return rows
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        groups[(r["experiment_accession"], r["file_format"])].append(r)
    out = []
    for grp in groups.values():
        if mode == "preferred":
            pref = [r for r in grp if r["preferred_default"]]
            out.extend(pref or grp)
        elif mode == "latest":
            out.append(max(grp, key=lambda r: r["date_created"]))
        else:
            raise ValueError(mode)
    return out


def download_one(sess: requests.Session, row: dict, outdir: str, verify_md5: bool) -> tuple[str, str]:
    subdir = row["experiment_accession"] or "_no_experiment"
    dest_dir = os.path.join(outdir, subdir)
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, os.path.basename(row["href"]))

    if os.path.exists(dest):
        if not verify_md5 or not row["md5sum"] or md5_of(dest) == row["md5sum"]:
            return dest, "exists"
        log(f"md5 mismatch for existing {dest}; re-downloading")

    tmp = dest + ".part"
    for attempt in range(1, 4):
        try:
            with sess.get(urljoin(BASE, row["href"]), stream=True, timeout=300) as r:
                r.raise_for_status()
                with open(tmp, "wb") as w:
                    for chunk in r.iter_content(1 << 20):
                        w.write(chunk)
            if verify_md5 and row["md5sum"] and md5_of(tmp) != row["md5sum"]:
                raise IOError("md5 mismatch")
            os.replace(tmp, dest)
            return dest, "downloaded"
        except Exception as e:  # noqa: BLE001
            log(f"{row['file_accession']} attempt {attempt} failed: {e}")
            time.sleep(2 * attempt)
    if os.path.exists(tmp):
        os.remove(tmp)
    return dest, "failed"


def step_extract(search_json: str, outdir: str = ".", formats: list[str] | None = None,
                 output_type: str = DEFAULT_OUTPUT_TYPE, dedup: str = "preferred",
                 download: bool = False, download_dir: str | None = None,
                 jobs: int = 4, verify_md5: bool = True,
                 targets: list[str] | None = None) -> str:
    formats = set(formats or ["bed"])
    download_dir = download_dir or os.path.join(outdir, "download")
    os.makedirs(outdir, exist_ok=True)

    with open(search_json) as f:
        records = json.load(f).get("@graph", [])

    rows, dropped = [], defaultdict(int)
    for rec in records:
        row = record_to_row(rec)
        if rec.get("status") != "released":
            dropped["not_released"] += 1; continue
        if row["output_type"] != output_type:
            dropped["output_type"] += 1; continue
        if row["file_format"] not in formats:
            dropped["file_format"] += 1; continue
        if not row["target"]:
            dropped["no_target"] += 1; continue
        if not row["href"]:
            dropped["no_href"] += 1; continue
        if targets and row["target"] not in targets:
            dropped["target_not_requested"] += 1; continue
        rows.append(row)

    n_before = len(rows)
    rows = dedup_rows(rows, dedup)
    dropped[f"dedup_{dedup}"] = n_before - len(rows)
    rows.sort(key=lambda r: (r["target"], r["experiment_accession"], r["file_accession"]))
    log(f"{len(records)} records -> {len(rows)} kept; dropped: {dict(dropped)}")

    for r in rows:
        sub = r["experiment_accession"] or "_no_experiment"
        r["local_path"] = os.path.abspath(os.path.join(download_dir, sub, os.path.basename(r["href"])))
        r["download_status"] = "exists" if os.path.exists(r["local_path"]) else "not_downloaded"

    if download:
        sess = session()
        todo = [r for r in rows]
        log(f"downloading {len(todo)} files to {download_dir} with {jobs} threads")
        done = 0
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(download_one, sess, r, download_dir, verify_md5): r for r in todo}
            for fut in as_completed(futs):
                r = futs[fut]
                path, status = fut.result()
                r["local_path"], r["download_status"] = os.path.abspath(path), status
                done += 1
                if done % 50 == 0 or done == len(todo):
                    log(f"  {done}/{len(todo)}")
        n_fail = sum(r["download_status"] == "failed" for r in rows)
        if n_fail:
            log(f"WARNING: {n_fail} downloads failed (see download_status in manifest); re-run to retry")

    manifest = os.path.join(outdir, "manifest.tsv")
    write_tsv(rows, manifest, MANIFEST_COLS)
    log(f"manifest -> {manifest}")
    return manifest


# ----------------------------------------------------------------------------
# step 3: stats
# ----------------------------------------------------------------------------
def step_stats(manifest: str, outdir: str = ".", compare_output_types: bool = False,
               assembly: str | None = None, rfa: str | None = "ENCODE4") -> dict:
    rows = read_tsv(manifest)
    by_tf: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_tf[r["target"]].append(r)

    tf_rows = []
    for tf, rs in sorted(by_tf.items()):
        bios = sorted({r["biosample"] for r in rs if r["biosample"]})
        tf_rows.append({
            "target": tf,
            "n_experiments": len({r["experiment_accession"] for r in rs}),
            "n_files": len(rs),
            "n_files_local": sum(r["download_status"] in ("exists", "downloaded") for r in rs),
            "n_biosamples": len(bios),
            "biosamples": ";".join(bios),
        })
    tf_tsv = os.path.join(outdir, "stats_per_target.tsv")
    write_tsv(tf_rows, tf_tsv, list(tf_rows[0].keys()) if tf_rows else ["target"])

    bios_count: dict[str, set] = defaultdict(set)
    for r in rows:
        bios_count[r["biosample"]].add(r["target"])
    bio_rows = [{"biosample": b, "n_targets": len(t)} for b, t in
                sorted(bios_count.items(), key=lambda kv: -len(kv[1]))]
    write_tsv(bio_rows, os.path.join(outdir, "stats_per_biosample.tsv"), ["biosample", "n_targets"])

    n_exp = [t["n_experiments"] for t in tf_rows]
    summary = {
        "n_files": len(rows),
        "n_experiments": len({r["experiment_accession"] for r in rows}),
        "n_targets": len(tf_rows),
        "n_biosamples": len(bio_rows),
        "targets_with_1_experiment": sum(x == 1 for x in n_exp),
        "targets_with_ge3_experiments": sum(x >= 3 for x in n_exp),
        "max_experiments_per_target": max(n_exp) if n_exp else 0,
    }

    # Like the old stats_tf_based.py: how many targets have conservative vs optimal IDR peaks
    if compare_output_types:
        if not assembly:
            assembly = rows[0]["assembly"] if rows else "mm10"
        types = ["conservative IDR thresholded peaks", "IDR thresholded peaks"]
        url = build_search_url(assembly, types, rfa, file_formats=["bed"],
                               fields=["target.label", "output_type"])
        log(f"GET {url}")
        r = session().get(url, timeout=600)
        r.raise_for_status()
        tf_to_ot: dict[str, set] = defaultdict(set)
        for rec in r.json().get("@graph", []):
            tf, ot = get_nested(rec, "target.label"), rec.get("output_type")
            if tf and ot:
                tf_to_ot[tf].add(ot)
        summary["coverage"] = {
            "all_targets_with_either": len(tf_to_ot),
            **{f"targets_with_{t}": sum(t in s for s in tf_to_ot.values()) for t in types},
            "only_optimal_IDR": sorted(tf for tf, s in tf_to_ot.items() if s == {types[1]}),
        }

    with open(os.path.join(outdir, "stats_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    for k, v in summary.items():
        if k != "coverage":
            log(f"{k}: {v}")
    if "coverage" in summary:
        c = summary["coverage"]
        log(f"coverage: {c['all_targets_with_either']} targets; "
            f"conservative={c['targets_with_conservative IDR thresholded peaks']}, "
            f"optimal={c['targets_with_IDR thresholded peaks']}, "
            f"only-optimal={len(c['only_optimal_IDR'])}")
    log(f"stats -> {outdir}/stats_*.tsv, stats_summary.json")
    return summary


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_search_args(sp):
        sp.add_argument("--assembly", default="mm10", help="mm10 | GRCh38 | ... (default mm10)")
        sp.add_argument("--rfa", default="ENCODE4", help="award.rfa filter; 'none' to disable")
        sp.add_argument("--output-type", default=DEFAULT_OUTPUT_TYPE)

    def add_extract_args(sp):
        sp.add_argument("--formats", nargs="+", default=["bed"],
                        help="file formats to keep (default: bed; bigBed holds the same peaks)")
        sp.add_argument("--dedup", choices=["preferred", "latest", "none"], default="preferred")
        sp.add_argument("--targets", nargs="*", help="restrict to these target labels")
        sp.add_argument("--download", action="store_true")
        sp.add_argument("--download-dir", help="default: <outdir>/download")
        sp.add_argument("-j", "--jobs", type=int, default=4)
        sp.add_argument("--no-md5", action="store_true", help="skip md5 verification")

    sp = sub.add_parser("search", help="step 1: query ENCODE, save raw JSON")
    sp.add_argument("--outdir", default=".")
    add_search_args(sp)

    sp = sub.add_parser("extract", help="step 2: filter -> manifest.tsv (+download)")
    sp.add_argument("search_json")
    sp.add_argument("--outdir", default=".")
    sp.add_argument("--output-type", default=DEFAULT_OUTPUT_TYPE)
    add_extract_args(sp)

    sp = sub.add_parser("stats", help="step 3: summary tables")
    sp.add_argument("manifest")
    sp.add_argument("--outdir", default=".")
    sp.add_argument("--compare-output-types", action="store_true",
                    help="also query ENCODE for conservative vs optimal IDR coverage")
    sp.add_argument("--assembly")
    sp.add_argument("--rfa", default="ENCODE4")

    sp = sub.add_parser("all", help="run search -> extract -> stats")
    sp.add_argument("--outdir", default=".")
    add_search_args(sp)
    add_extract_args(sp)
    sp.add_argument("--compare-output-types", action="store_true")

    a = p.parse_args(argv)
    rfa = None if getattr(a, "rfa", None) in (None, "none", "") else a.rfa

    if a.cmd == "search":
        step_search(a.assembly, a.outdir, rfa, [a.output_type])
    elif a.cmd == "extract":
        step_extract(a.search_json, a.outdir, a.formats, a.output_type, a.dedup,
                     a.download, a.download_dir, a.jobs, not a.no_md5, a.targets)
    elif a.cmd == "stats":
        step_stats(a.manifest, a.outdir, a.compare_output_types, a.assembly, rfa)
    elif a.cmd == "all":
        js = step_search(a.assembly, a.outdir, rfa, [a.output_type])
        man = step_extract(js, a.outdir, a.formats, a.output_type, a.dedup,
                           a.download, a.download_dir, a.jobs, not a.no_md5, a.targets)
        step_stats(man, a.outdir, a.compare_output_types, a.assembly, rfa)


if __name__ == "__main__":
    main()
