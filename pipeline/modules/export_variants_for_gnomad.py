#!/usr/bin/env python
import pandas as pd
import argparse
import os

# Columnas que produce summarize_*_vcf (sin header)
COLS = ["CHROM","POS","REF","ALT",
        "SAMPLE_COUNT","CLINVAR_MATCH","CLINVAR_INFO","SOURCE"]

def norm_chr(c):
    c = str(c)
    if c.lower().startswith("chr"):
        c = c[3:]
    if c in ["23","X","x"]:
        return "X"
    if c in ["24","Y","y"]:
        return "Y"
    if c.upper() in ["MT","M","MITO","MTE","MTE?"]:
        return "MT"
    return c

def main(args):
    som_path = args.somatic_tsv or "data/processed/somatic_summary.tsv"
    ger_path = args.germline_tsv or "data/processed/germline_summary.tsv"
    out_path = args.out_tsv or "data/processed/variants_for_gnomad_query.tsv"

    # Leer TSVs SIN header (como salen de bcftools query)
    som  = pd.read_csv(som_path,  sep="\t", header=None, names=COLS, dtype=str)
    germ = pd.read_csv(ger_path,  sep="\t", header=None, names=COLS, dtype=str)

    both = pd.concat([som, germ], ignore_index=True)

    both["CHR_NORM"] = both["CHROM"].apply(norm_chr)

    # variant_id que espera gnomAD (GRCh38)
    both["variant_id"] = (
        both["CHR_NORM"].astype(str) + "-" +
        both["POS"].astype(str) + "-" +
        both["REF"].astype(str) + "-" +
        both["ALT"].astype(str)
    )

    unique_for_gnomad = (
        both[["CHR_NORM","POS","REF","ALT","variant_id"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    unique_for_gnomad.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] Variantes únicas listas para gnomAD -> {out_path}")
    print(unique_for_gnomad.head())

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--somatic-tsv",   default=None)
    p.add_argument("--germline-tsv",  default=None)
    p.add_argument("--out-tsv",       default=None)
    main(p.parse_args())
