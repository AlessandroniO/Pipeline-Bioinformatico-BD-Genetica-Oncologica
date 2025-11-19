#!/usr/bin/env python3
import sys
import os
import pandas as pd

INFO_FIELDS = ["SAMPLE_COUNT", "CLINVAR_MATCH", "SOURCE", "CLINVAR_INFO"]

CONTIGS = (["chr" + str(i) for i in range(1, 23)] + ["chr23", "chrY", "chrMT"])

def chrom_rank(ch):
    order = {c: i for i, c in enumerate(CONTIGS)}
    return order.get(str(ch), 10**9)

def norm_chr_to_vcf(ch):
    if pd.isna(ch):
        return None
    s = str(ch).strip()
    if s == "":
        return None
    sl = s.lower()
    if sl.startswith("chr"):
        s = s
    elif sl.startswith("nc_"):
        core = s.replace("NC_", "").split(".")[0]
        core = core.lstrip("0") or "0"
        s = "chr" + core
    elif s.isdigit():
        s = "chr" + s
    elif s in ["X", "x", "23"]:
        s = "chr23"
    elif s in ["Y", "y", "24"]:
        s = "chrY"
    elif s.upper() in ["MT", "M", "MTDNA"]:
        s = "chrMT"
    else:
        if not s.lower().startswith("chr"):
            s = "chr" + s
    return s

def pick_coords(df):
    df = df.copy()
    cands = [
        ("CHROM_norm", "POS_norm", "REF_from_hgvs", "ALT_from_hgvs"),
        ("#CHROM",     "POS",      "REF",           "ALT"),
        ("CHROM_from_hgvs", "POS_from_hgvs", "REF_from_hgvs", "ALT_from_hgvs"),
    ]
    df["_CHROM"] = pd.NA
    df["_POS"]   = pd.NA
    df["_REF"]   = pd.NA
    df["_ALT"]   = pd.NA
    for c_chrom, c_pos, c_ref, c_alt in cands:
        have = all(col in df.columns for col in [c_chrom, c_pos, c_ref, c_alt])
        if not have:
            continue
        mask = (
            df["_CHROM"].isna()
            & df[c_chrom].notna()
            & df[c_pos].notna()
            & df[c_ref].notna()
            & df[c_alt].notna()
        )
        df.loc[mask, "_CHROM"] = df.loc[mask, c_chrom]
        df.loc[mask, "_POS"]   = df.loc[mask, c_pos]
        df.loc[mask, "_REF"]   = df.loc[mask, c_ref].astype(str)
        df.loc[mask, "_ALT"]   = df.loc[mask, c_alt].astype(str)

    df["_CHROM"] = df["_CHROM"].apply(norm_chr_to_vcf)
    df["_POS"]   = pd.to_numeric(df["_POS"], errors="coerce").astype("Int64")
    df["_REF"]   = df["_REF"].astype("string")
    df["_ALT"]   = df["_ALT"].astype("string")

    df_out = df.dropna(subset=["_CHROM", "_POS", "_REF", "_ALT"]).copy()
    df_out.rename(columns={"_CHROM": "CHROM", "_POS": "POS", "_REF": "REF", "_ALT": "ALT"}, inplace=True)
    return df_out

def write_vcf(df, label, out_vcf):
    os.makedirs(os.path.dirname(out_vcf), exist_ok=True)
    header = []
    header.append("##fileformat=VCFv4.2")
    header.append("##FILTER=<ID=PASS,Description=\"All filters passed\">")
    header.append("##reference=GRCh38")
    header.append("##INFO=<ID=SAMPLE_COUNT,Number=1,Type=Integer,Description=\"How many patient-level rows carried this variant in this cohort\">")
    header.append("##INFO=<ID=CLINVAR_MATCH,Number=1,Type=Integer,Description=\"1 if variant matched ClinVar entry in join_with_vcf step\">")
    header.append("##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"SOMATIC or GERMINAL cohort label\">")
    header.append("##INFO=<ID=CLINVAR_INFO,Number=1,Type=String,Description=\"ClinVar INFO_STRING subset for this variant (first hit)\">")
    for ctg in CONTIGS:
        header.append(f"##contig=<ID={ctg}>")
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    df = df.copy()
    df["chrom_order"] = df["CHROM"].apply(chrom_rank)
    df = df.sort_values(["chrom_order", "POS", "REF", "ALT"], kind="mergesort")
    df = df.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"])

    for c in INFO_FIELDS:
        if c not in df.columns:
            df[c] = pd.NA

    def build_info(row):
        parts = []
        if pd.notna(row.get("SAMPLE_COUNT", pd.NA)):
            parts.append(
                f"SAMPLE_COUNT={int(row['SAMPLE_COUNT'])}"
                if str(row["SAMPLE_COUNT"]).isdigit()
                else f"SAMPLE_COUNT={row['SAMPLE_COUNT']}"
            )
        if pd.notna(row.get("CLINVAR_MATCH", pd.NA)):
            parts.append(f"CLINVAR_MATCH={row['CLINVAR_MATCH']}")
        if pd.notna(row.get("SOURCE", pd.NA)):
            parts.append(f"SOURCE={row['SOURCE']}")
        if pd.notna(row.get("CLINVAR_INFO", pd.NA)) and str(row["CLINVAR_INFO"]).strip() != "":
            v = str(row["CLINVAR_INFO"]).replace("\t", " ").replace("\n", " ")
            parts.append(f"CLINVAR_INFO={v}")
        return ";".join(parts) if parts else "."

    df["INFO"] = df.apply(build_info, axis=1)

    with open(out_vcf, "w", encoding="utf-8") as f:
        for line in header:
            f.write(line + "\n")
        for _, r in df.iterrows():
            f.write(f"{r['CHROM']}\t{int(r['POS'])}\t.\t{r['REF']}\t{r['ALT']}\t.\tPASS\t{r['INFO']}\n")

    print(f"[OK] VCF escrito en {out_vcf}")
    print(f"Total variantes únicas: {len(df)}")

def ensure_series(df, colname):
    """Si df[colname] es DataFrame (duplicados), toma la 1ª no-nula por fila."""
    obj = df[colname]
    if isinstance(obj, pd.DataFrame):
        s = obj.bfill(axis=1).iloc[:, 0]
    else:
        s = obj
    return s.astype("string")

def main(in_csv, label, out_vcf):
    df = pd.read_csv(in_csv, low_memory=False)

    # 1) Coordenadas desde el CSV
    vcf_df = pick_coords(df)

    # 2) Normaliza nombres y elimina duplicados EXACTOS de columnas
    vcf_df.columns = vcf_df.columns.map(str).str.strip()
    vcf_df = vcf_df.loc[:, ~vcf_df.columns.duplicated()].copy()

    if vcf_df.empty:
        raise RuntimeError(
            "No hay ninguna fila con coordenadas completas para VCF. "
            "Revisá que existan columnas CHROM_norm/POS_norm/REF_from_hgvs/ALT_from_hgvs "
            "o bien #CHROM/POS/REF/ALT en el CSV de entrada."
        )

    # 3) Asegura que REF/ALT sean Series aun si hay duplicados implícitos
    #    (p.ej., merges que dejaron varias columnas llamadas 'REF'/'ALT')
    vcf_df["REF"] = ensure_series(vcf_df, "REF")
    vcf_df["ALT"] = ensure_series(vcf_df, "ALT")

    # 4) Filtro SNV (1 base A/C/G/T)
    before = len(vcf_df)
    vcf_df = vcf_df[
        vcf_df["REF"].str.fullmatch(r"[ACGT]", na=False) &
        vcf_df["ALT"].str.fullmatch(r"[ACGT]", na=False)
    ].copy()
    print(f"[INFO] Variantes con coords: {before} | SNVs válidas para VCF: {len(vcf_df)}")

    # 5) Copia campos INFO si están en el CSV original
    for f in INFO_FIELDS:
        if f in df.columns and f not in vcf_df.columns:
            vcf_df[f] = df[f]

    write_vcf(vcf_df, label, out_vcf)

    # ----------------------------------------------------------------------

    if vcf_df.empty:
        raise RuntimeError(
            "No hay ninguna fila con coordenadas completas para VCF. "
            "Revisá que existan columnas CHROM_norm/POS_norm/REF_from_hgvs/ALT_from_hgvs "
            "o bien #CHROM/POS/REF/ALT en el CSV de entrada."
        )

    for f in INFO_FIELDS:
        if f in df.columns and f not in vcf_df.columns:
            vcf_df[f] = df[f]

    write_vcf(vcf_df, label, out_vcf)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python csv_to_vcf.py <input_csv> <SOMATIC|GERMINAL> <output_vcf>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
