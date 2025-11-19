#!/usr/bin/env python3
import os, sys, pandas as pd, unicodedata

SUMMARY = "data/processed/summary_somatic_tumor_gene.csv"
CIE10   = "Lista_tumores_CIE-10.csv"
OUT     = "data/processed/gene_tissue_map.csv"

def strip_accents(s: str) -> str:
    return "".join(c for c in unicodedata.normalize("NFKD", s) if not unicodedata.combining(c))

def normcols(df):
    # minúsculas, sin espacios y sin acentos
    df.columns = [strip_accents(c.strip().lower()) for c in df.columns]
    return df

def pick_col(df, candidates, required_name):
    env_override = os.environ.get(required_name.upper().replace(" ", "_"))  # e.g. TUMOR_GROUP_(EN_EL_MAPEO)
    if env_override:
        key = strip_accents(env_override.strip().lower())
        if key in df.columns:
            return key
        else:
            raise RuntimeError(f"El override de entorno '{env_override}' no existe en columnas {list(df.columns)}")
    for c in candidates:
        if c in df.columns:
            return c
    raise RuntimeError(
        f"No se encontró la columna requerida '{required_name}'. "
        f"Probé con alias: {candidates} en {list(df.columns)}"
    )

def key_norm(s):
    return str(s).strip().lower()

def main():
    # --- Leer summary somático (ya generado por tu pipeline)
    df = pd.read_csv(SUMMARY)
    df = normcols(df)
    need = {"tumor_group", "gene", "n_filas"}
    if not need.issubset(df.columns):
        raise RuntimeError(f"En {SUMMARY} faltan columnas {need - set(df.columns)}")

    # --- Leer tu mapeo CIE-10 → tejido
    mapdf = pd.read_csv(CIE10)
    mapdf = normcols(mapdf)

    # Alias más amplios (incluye ddescripcion / descripcion)
    tumor_candidates = [
        "tumor_group", "grupo_tumor", "grupo", "tumor",
        "diagnostico", "diagnostico_principal", "dx",
        "cie10_desc", "descripcion", "ddescripcion", "desc", "nombre", "nombre_tumor",
        "cie10", "cie10_descripcion", "detalle"
    ]
    tissue_candidates = [
        "tissue", "tejido", "organ", "organo", "sistema", "system", "site", "sitio",
        "ubicacion", "localizacion"
    ]

    # OJO: los nombres para overrides por entorno son:
    #   TUMOR_COL  -> fuerza cuál usar como tumor_group en el mapeo
    #   TISSUE_COL -> fuerza cuál usar como tissue en el mapeo
    tumor_env = os.environ.get("TUMOR_COL")
    tissue_env = os.environ.get("TISSUE_COL")

    col_tumor_map = strip_accents(tumor_env.strip().lower()) if tumor_env else pick_col(mapdf, tumor_candidates, "tumor_group (en el mapeo)")
    if tissue_env:
        col_tissue = strip_accents(tissue_env.strip().lower())
        if col_tissue not in mapdf.columns:
            raise RuntimeError(f"El override TISSUE_COL='{tissue_env}' no existe en columnas {list(mapdf.columns)}")
    else:
        # si no hay tissue, usamos el mismo tumor_group como proxy
        try:
            col_tissue = pick_col(mapdf, tissue_candidates, "tissue/tejido (en el mapeo)")
        except RuntimeError:
            col_tissue = col_tumor_map

    print(f"[INFO] Columna de mapeo (tumor_group): {col_tumor_map}")
    print(f"[INFO] Columna de mapeo (tissue)     : {col_tissue}")

    # Normalizar claves para el join
    df["tumor_group_key"]    = df["tumor_group"].apply(key_norm)
    mapdf["tumor_group_key"] = mapdf[col_tumor_map].apply(key_norm)

    # Join (left) de summary contra el diccionario
    j = pd.merge(
        df,
        mapdf[["tumor_group_key", col_tissue]].rename(columns={col_tissue: "tissue"}),
        on="tumor_group_key",
        how="left"
    )
    j["tissue"] = j["tissue"].fillna("UNKNOWN")

    # Agregado final por gen y tejido
    agg = (
        j.groupby(["gene", "tissue"], as_index=False)
         .agg(
            n_var_cohort=("n_filas", "sum"),
            tumor_groups=("tumor_group", lambda s: sorted(set(s)))
         )
    )
    agg["n_tumor_groups"] = agg["tumor_groups"].apply(len)
    agg = agg.sort_values(["gene", "n_var_cohort"], ascending=[True, False])

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    agg.to_csv(OUT, index=False)
    print(f"[OK] {OUT} escrito. Filas: {len(agg)}")

if __name__ == "__main__":
    main()
