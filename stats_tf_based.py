import requests
from collections import defaultdict

url = (
    "https://www.encodeproject.org/search/?type=File"
    "&assay_title=TF+ChIP-seq"
    "&assembly=GRCh38"
    "&status=released"
    "&file_format=bed"
    "&award.rfa=ENCODE4"
    "&field=output_type"
    "&field=accession"
    "&field=target.label"
    "&limit=all"
    "&format=json"
)
headers = {"accept": "application/json"}

r = requests.get(url, headers=headers)
r.raise_for_status()
files = r.json()["@graph"]

tf_to_outputs = defaultdict(set)

for f in files:
    tf = (f.get("target", {}) or {}).get("label")
    ot = f.get("output_type")
    if not tf or not ot:
        continue
    tf_to_outputs[tf].add(ot)

TFs = set(tf_to_outputs.keys())

want_conservative = "conservative IDR thresholded peaks"
want_idr = "IDR thresholded peaks"

tfs_conservative = {tf for tf, outs in tf_to_outputs.items() if want_conservative in outs}
tfs_idr = {tf for tf, outs in tf_to_outputs.items() if want_idr in outs}

print(f"Total TF labels: {len(TFs)}")
print(f"TFs with '{want_conservative}': {len(tfs_conservative)}")
print(f"TFs with '{want_idr}': {len(tfs_idr)}")