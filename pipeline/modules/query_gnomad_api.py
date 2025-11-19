#!/usr/bin/env python
import pandas as pd
import requests
import time
import argparse
import os

ENDPOINT = "https://gnomad.broadinstitute.org/api"
GRAPHQL_QUERY = """
query ($variantId: String!) {
  variant(variantId: $variantId, dataset: gnomad_r4) {
    variantId
    genome { ac an af }
  }
}
"""

def main(args):
    in_path  = args.in_tsv  or "data/processed/variants_for_gnomad_query.tsv"
    out_path = args.out_tsv or "data/processed/gnomad_results.tsv"

    dfq = pd.read_csv(in_path, sep="\t", dtype=str)
    results = []

    for vid in dfq["variant_id"]:
        payload = {"query": GRAPHQL_QUERY, "variables": {"variantId": vid}}
        try:
            r = requests.post(ENDPOINT, json=payload, timeout=10)
            r.raise_for_status()
            js = r.json()
            v = (js.get("data") or {}).get("variant")
            if not v:
                results.append({"variant_id": vid, "ac": None, "an": None, "af": None})
            else:
                g = v.get("genome") or {}
                results.append({
                    "variant_id": v.get("variantId", vid),
                    "ac": g.get("ac"),
                    "an": g.get("an"),
                    "af": g.get("af"),
                })
        except Exception:
            results.append({"variant_id": vid, "ac": None, "an": None, "af": None})
        time.sleep(0.2)

    out_df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] Frecuencias gnomAD guardadas en {out_path}")
    print(out_df.head())

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--in-tsv",  default=None)
    p.add_argument("--out-tsv", default=None)
    main(p.parse_args())
