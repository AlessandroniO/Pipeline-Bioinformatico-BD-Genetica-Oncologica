import os
import re
import pandas as pd
import numpy as np


# Base del proyecto: por defecto /app (dentro del contenedor).
# En Windows local podrías exportar BASE_DIR, pero no es necesario si corrés en Docker.
BASE = os.environ.get("BASE_DIR", "/app")

# Rutas portables (sin backslashes)
SOMATIC_CSV   = os.path.join(BASE, "data", "intermediate", "variantes_somaticas_API_anotadas.csv")
GERMLINE_CSV  = os.path.join(BASE, "data", "intermediate", "variantes_germinales_API_anotadas.csv")

OUT_TOP_SOM   = os.path.join(BASE, "data", "processed", "summary_top_genes_somatic.csv")
OUT_TOP_GER   = os.path.join(BASE, "data", "processed", "summary_top_genes_germline.csv")
OUT_ACC_SOM   = os.path.join(BASE, "data", "processed", "summary_accionables_somatic.csv")
OUT_TUMOR_MAP = os.path.join(BASE, "data", "processed", "summary_somatic_tumor_gene.csv")

os.makedirs(os.path.join(BASE, r"data\processed"), exist_ok=True)

GENE_COL_CANDIDATES = [
    "gene", "GENE", "Gene", "gen", "symbol", "Symbol"
]

INFO_COL_CANDIDATES = [
    "INFO_STRING", "info_string", "INFO", "ClinVar_INFO", "CLINVAR_INFO"
]

def first_existing_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def infer_gene_from_variante(series_variante: pd.Series) -> pd.Series:
    """
    Intenta inferir el símbolo de gen buscando el PRIMER texto entre paréntesis.
    Ej:
      "c.2098A>G (SF3B1) (p.Lys700Glu)"  -> SF3B1
      "PIK3CA (c.1633G>A) (p.Glu545Lys)" -> PIK3CA
    Esta versión NO asume que el gen está al principio.
    """
    pat = re.compile(r"\(([A-Za-z0-9._-]+)\)")
    inferred = []
    for val in series_variante.astype(str):
        m = pat.search(val)
        if m:
            inferred.append(m.group(1))
        else:
            m2 = re.match(r"^([A-Za-z0-9._-]+)\s*\(", val)
            inferred.append(m2.group(1) if m2 else np.nan)
    return pd.Series(inferred, index=series_variante.index, dtype="object")

def clean_gene_symbol(g):
    if pd.isna(g):
        return np.nan
    g = str(g).strip()
    if g == "":
        return np.nan
    if g.lower().startswith("chr") and ":" in g:
        return np.nan
    return g

def mark_actionable(info_str):
    if pd.isna(info_str):
        return False
    text = str(info_str).lower()
    keywords = [
        "pathogenic",
        "likely_pathogenic",
        "drug",
        "therapy",
        "therapeutic",
        "sensitivity",
        "response",
        "target",
    ]
    return any(k in text for k in keywords)

def build_top_genes_table(df_in: pd.DataFrame, cohort_label: str) -> pd.DataFrame:
    gene_col = first_existing_col(df_in, GENE_COL_CANDIDATES)
    if gene_col is None:
        if "variante" in df_in.columns:
            df_in["_gene_tmp"] = infer_gene_from_variante(df_in["variante"])
            gene_col = "_gene_tmp"
        else:
            df_in["_gene_tmp"] = np.nan
            gene_col = "_gene_tmp"

    df_in["_gene_clean"] = df_in[gene_col].apply(clean_gene_symbol)

    counts = (
        df_in
        .groupby("_gene_clean", dropna=True)
        .size()
        .reset_index(name="n_filas")
        .sort_values("n_filas", ascending=False)
        .head(30)
        .reset_index(drop=True)
    )

    counts.rename(columns={"_gene_clean": "gene"}, inplace=True)
    counts["cohorte"] = cohort_label.upper()
    return counts

def build_actionable_table(df_somatic: pd.DataFrame) -> pd.DataFrame:
    gene_col = first_existing_col(df_somatic, GENE_COL_CANDIDATES)
    if gene_col is None:
        if "variante" in df_somatic.columns:
            df_somatic["_gene_tmp"] = infer_gene_from_variante(df_somatic["variante"])
            gene_col = "_gene_tmp"
        else:
            df_somatic["_gene_tmp"] = np.nan
            gene_col = "_gene_tmp"

    info_col = first_existing_col(df_somatic, INFO_COL_CANDIDATES)
    if info_col is None:
        return pd.DataFrame(
            columns=["gene","CHROM","POS","REF","ALT","INFO_STRING","es_accionable"]
        )

    chrom_col = None
    pos_col   = None
    ref_col   = None
    alt_col   = None

    for cand in ["CHROM_norm", "CHROM_from_hgvs", "CHROM", "#CHROM", "chr", "chrom"]:
        if cand in df_somatic.columns:
            chrom_col = cand
            break
    for cand in ["POS_norm", "POS_from_hgvs", "POS", "pos"]:
        if cand in df_somatic.columns:
            pos_col = cand
            break
    for cand in ["REF_from_hgvs", "REF", "ref"]:
        if cand in df_somatic.columns:
            ref_col = cand
            break
    for cand in ["ALT_from_hgvs", "ALT", "alt"]:
        if cand in df_somatic.columns:
            alt_col = cand
            break

    subcols = {
        "gene": df_somatic[gene_col].apply(clean_gene_symbol),
        "INFO_STRING": df_somatic[info_col].astype(str),
        "es_accionable": df_somatic[info_col].apply(mark_actionable),
    }

    if chrom_col: subcols["CHROM"] = df_somatic[chrom_col]
    if pos_col:   subcols["POS"]   = df_somatic[pos_col]
    if ref_col:   subcols["REF"]   = df_somatic[ref_col]
    if alt_col:   subcols["ALT"]   = df_somatic[alt_col]

    out = pd.DataFrame(subcols)
    out = out[out["es_accionable"] == True].copy()

    sort_cols = [c for c in ["gene","CHROM","POS"] if c in out.columns]
    out = out.sort_values(by=sort_cols, ascending=[True]*len(sort_cols), na_position="last")
    out.reset_index(drop=True, inplace=True)

    return out

def build_tumor_gene_table_light(df_somatic: pd.DataFrame) -> pd.DataFrame:
    tumor_col_candidates = [
        "clasificacion_final",
        "razon_clasificacion",
        "tipo_muestra_norm",
        "tipo_muestra_normalizada"
    ]

    tumor_group = None
    for c in tumor_col_candidates:
        if c in df_somatic.columns:
            if tumor_group is None:
                tumor_group = df_somatic[c].astype(str)
            else:
                tumor_group = tumor_group.fillna(df_somatic[c].astype(str))

    if tumor_group is None:
        tumor_group = pd.Series([np.nan]*len(df_somatic), index=df_somatic.index)

    tumor_blocklist = [
        "nan","none","sangre","blood","plasma",
        "líquido pleural","liquido pleural",
        "tumor en parafina (ffpe)","tumor en parafina","ffpe"
    ]
    tumor_group = tumor_group.str.strip().str.lower()
    tumor_group = tumor_group.where(~tumor_group.isin(tumor_blocklist))

    df_somatic = df_somatic.copy()
    df_somatic["tumor_group"] = tumor_group

    if "variante" in df_somatic.columns:
        df_somatic["gene_tmp"] = infer_gene_from_variante(df_somatic["variante"])
    else:
        df_somatic["gene_tmp"] = np.nan

    df_somatic["gene_tmp"] = df_somatic["gene_tmp"].apply(clean_gene_symbol)

    df_valid = df_somatic.dropna(subset=["tumor_group","gene_tmp"]).copy()
    df_valid = df_valid[
        (df_valid["tumor_group"] != "") &
        (df_valid["gene_tmp"] != "")
    ].copy()

    if df_valid.empty:
        return pd.DataFrame(columns=["tumor_group","gene","n_filas","rank_in_tumor"])

    counts = (
        df_valid
        .groupby(["tumor_group","gene_tmp"])
        .size()
        .reset_index(name="n_filas")
    )

    counts["rank_in_tumor"] = counts.groupby("tumor_group")["n_filas"] \
                                    .rank(method="first", ascending=False)

    top_counts = counts[counts["rank_in_tumor"] <= 10].copy()
    top_counts = top_counts.sort_values(
        ["tumor_group","n_filas"],
        ascending=[True, False]
    ).reset_index(drop=True)

    top_counts = top_counts.rename(columns={"gene_tmp": "gene"})
    return top_counts[["tumor_group","gene","n_filas","rank_in_tumor"]]

def main():
    print(">> Leyendo tablas de entrada...")
    df_som = pd.read_csv(SOMATIC_CSV, low_memory=False)
    df_ger = pd.read_csv(GERMLINE_CSV, low_memory=False)

    print(">> Generando Paso A (top genes somática/germinal)...")
    top_som = build_top_genes_table(df_som.copy(), cohort_label="somatic")
    top_ger = build_top_genes_table(df_ger.copy(), cohort_label="germline")
    top_som.to_csv(OUT_TOP_SOM, index=False, encoding="utf-8")
    top_ger.to_csv(OUT_TOP_GER, index=False, encoding="utf-8")
    print(f"   Guardado A1: {OUT_TOP_SOM}")
    print(f"   Guardado A2: {OUT_TOP_GER}")

    print(">> Generando Paso B (variantes somáticas accionables)...")
    accionables = build_actionable_table(df_som.copy())
    accionables.to_csv(OUT_ACC_SOM, index=False, encoding="utf-8")
    print(f"   Guardado B : {OUT_ACC_SOM} (filas accionables = {len(accionables)})")

    print(">> Generando Paso C (mapa Tumor ↔ Gen, versión liviana)...")
    tumor_map_light = build_tumor_gene_table_light(df_som.copy())
    tumor_map_light.to_csv(OUT_TUMOR_MAP, index=False, encoding="utf-8")
    print(f"   Guardado C : {OUT_TUMOR_MAP} (filas = {len(tumor_map_light)})")

    print("\n=== RESUMEN LISTO ===")
    print("Tabla / Figura 1:")
    print(f"  {OUT_TOP_SOM}")
    print(f"  {OUT_TOP_GER}")
    print("Tabla / Figura 2:")
    print(f"  {OUT_ACC_SOM}")
    print("Tabla / Figura 3:")
    print(f"  {OUT_TUMOR_MAP}")

if __name__ == "__main__":
    main()

