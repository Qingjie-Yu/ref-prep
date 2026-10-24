#!/usr/bin/env python3
import json, requests, argparse, os, csv
from urllib.parse import urljoin

BASE = "https://www.encodeproject.org"
HEADERS = {"accept": "application/json"}

def get_json(path_or_url: str):
    # Accept either "/files/ENCFF.../" or a full URL
    url = path_or_url
    if url.startswith("/"):
        url = urljoin(BASE, url)
    if "?" not in url:
        url = url + "?format=json&frame=object"
    r = requests.get(url, headers=HEADERS)
    r.raise_for_status()
    return r.json()

def main():
    p = argparse.ArgumentParser()
    p.add_argument("file_search_json", help="JSON from ENCODE search/?type=File ...")
    p.add_argument("--download", action="store_true")
    p.add_argument("--out", default="./download")
    p.add_argument("--csv", default="encode_files.csv")
    p.add_argument("--formats", nargs="*", default=["bed", "bigBed"],
                   help="File formats to keep (default: bed bigBed)")
    args = p.parse_args()

    with open(args.file_search_json) as f:
        items = json.load(f).get("@graph", [])

    want_formats = set(args.formats)
    rows = []

    os.makedirs(args.out, exist_ok=True)

    for obj in items:
        # Some searches return partial objects; if href/output_type/etc missing, fetch full file object
        fj = obj
        if fj.get("href") is None or fj.get("file_format") is None or fj.get("output_type") is None:
            fid = fj.get("@id")
            if not fid:
                continue
            try:
                fj = get_json(fid)
            except Exception:
                continue

        if fj.get("status") != "released":
            continue
        if fj.get("output_type") != "conservative IDR thresholded peaks":
            continue
        if fj.get("file_format") not in want_formats:
            continue

        file_acc = fj.get("accession")
        href = fj.get("href")
        fmt = fj.get("file_format")
        target = (fj.get("target") or {}).get("label", "")

        # dataset is usually the experiment @id like "/experiments/ENCSR....../"
        dataset = fj.get("dataset") or ""
        exp_acc = ""
        if isinstance(dataset, str) and dataset.startswith("/experiments/"):
            exp_acc = dataset.strip("/").split("/")[1]

        rows.append({
            "experiment_accession": exp_acc,
            "target": target,
            "file_accession": file_acc,
            "href": href,
            "file_format": fmt,
            "output_type": fj.get("output_type"),
            "status": fj.get("status"),
        })

        if args.download and href:
            # put downloads under out/ENCSR.../ (if we have exp), else out/_no_experiment/
            subdir = exp_acc if exp_acc else "_no_experiment"
            out_dir = os.path.join(args.out, subdir)
            os.makedirs(out_dir, exist_ok=True)
            out_path = os.path.join(out_dir, os.path.basename(href))

            if not os.path.exists(out_path):
                with requests.get(urljoin(BASE, href), stream=True, headers=HEADERS) as r:
                    r.raise_for_status()
                    with open(out_path, "wb") as w:
                        for chunk in r.iter_content(1024 * 1024):
                            w.write(chunk)

    csv_path = os.path.join(args.out, args.csv)
    fieldnames = ["experiment_accession", "target", "file_accession", "href",
                  "file_format", "output_type", "status"]

    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {csv_path}")

if __name__ == "__main__":
    main()