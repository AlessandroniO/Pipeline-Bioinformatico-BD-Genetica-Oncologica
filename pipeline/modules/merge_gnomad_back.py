#!/usr/bin/env python
import pandas as pd
import argparse
import os

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
    som_in   = args.somatic_tsv or "data/processed/somatic_summary.tsv"
    ger_in   = args.germline_tsv or "data/processed/germline_summary.tsv"
    freqs_in = args.freqs_tsv   or "data/processed/gnomad_results.tsv"
    som_out  = args.som_out     or "data/processed/somatic_with_gnomad.tsv"
    ger_out  = args.germ_out    or "data/processed/germline_with_gnomad.tsv"

    # Si alguno no existe o está vacío, fallar explícito (no silencioso)
    for p in [som_in, ger_in, freqs_in]:
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            raise RuntimeError(f"Entrada inexistente o vacía: {p}")

    som  = pd.read_csv(som_in, sep="\t", header=None, names=COLS, dtype=str)
    germ = pd.read_csv(ger_in, sep="\t", header=None, names=COLS, dtype=str)
    freqs = pd.read_csv(freqs_in, sep="\t", dtype=str)

    for df in (som, germ):
        df["CHR_NORM"] = df["CHROM"].apply(norm_chr)
        df["variant_id"] = (
            df["CHR_NORM"].astype(str) + "-" +
            df["POS"].astype(str) + "-" +
            df["REF"].astype(str) + "-" +
            df["ALT"].astype(str)
        )

    som_annot  = som.merge(freqs, on="variant_id", how="left")
    germ_annot = germ.merge(freqs, on="variant_id", how="left")

    os.makedirs(os.path.dirname(som_out), exist_ok=True)
    som_annot.to_csv(som_out, sep="\t", index=False)
    germ_annot.to_csv(ger_out, sep="\t", index=False)

    print(f"[OK] somatic_with_gnomad.tsv -> {som_out}")
    print(f"[OK] germline_with_gnomad.tsv -> {ger_out}")
    print("Ejemplo somatic_with_gnomad:")
    print(som_annot.head())

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--somatic-tsv",   default=None)
    p.add_argument("--germline-tsv",  default=None)
    p.add_argument("--freqs-tsv",     default=None)
    p.add_argument("--som-out",       default=None)
    p.add_argument("--germ-out",      default=None)
    main(p.parse_args())
