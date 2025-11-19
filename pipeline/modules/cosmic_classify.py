# cosmic_classify.py
# Requiere: pandas, duckdb, pyarrow, numpy, regex
#
# Objetivo:
#   - Leer variantes_extraidas.csv (todas las filas paciente-variantes)
#   - Cruzar con COSMIC (archivo TSV enorme, local)
#   - Inferir evidencia somática por COSMIC
#   - Combinar eso con tipo_de_muestra (sangre / tumor / etc.)
#   - Devolver un CSV enriquecido:
#       es_reportada_somatica        (True/False/NaN a nivel id)
#       clasificacion_final          ('somática', 'somática (solo por origen tumoral)',
#                                     'germinal', 'indeterminado')
#       razon_clasificacion          (explicación corta)

import sys
import re
import numpy as np
import pandas as pd
import duckdb

# ---------- helpers internos ----------

pat_gene  = re.compile(r"\(([^)]+)\)")
pat_c     = re.compile(r"(c\.[^ )]+)")
pat_p     = re.compile(r"\((p\.[^)]+)\)")

def extract_first(pattern, x):
    if pd.isna(x):
        return np.nan
    m = pattern.search(str(x))
    return m.group(1) if m else np.nan

# mapa 3 letras -> 1 letra
aa3to1 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C",
    "Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I",
    "Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P",
    "Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
    "Ter":"*"
}
valid_pat = re.compile(r"^[A-Za-z]{3}[0-9]+[A-Za-z=*]{0,3}$")

def aa3_to_aa1(s):
    if pd.isna(s):
        return np.nan
    s = str(s)
    if not valid_pat.match(s):
        return np.nan
    three = s[:3]
    rest  = s[3:]
    one   = aa3to1.get(three, three)  # fallback: dejar igual si no conocemos
    return f"{one}{rest}"

def reduce_flag(series):
    """
    series es es_reportada_somatica por un mismo id (misma variante en ese informe).
    Regla:
    - Si hay al menos un True -> True
    - Si TODO es NaN          -> NaN
    - Sino                    -> False
    """
    if series.any(skipna=True):
        return True
    if series.isna().all():
        return np.nan
    return False

def normalizar_tipo_muestra(tm_raw: str):
    """
    Devuelve una etiqueta canónica:
      'sangre'  -> "sangre"
      'tumor', 'tumor en parafina (ffpe)', 'líquido pleural' -> "tumor"
      else -> "desconocida"
    Trabajamos todo lowercase para robustez.
    """
    if tm_raw is None:
        return "desconocida"
    t = str(tm_raw).strip().lower()

    # corregimos encoding raro de "líquido"
    t = t.replace("lÃ­quido", "líquido")

    # mapeo explícito
    if t in ["sangre"]:
        return "sangre"
    if t in [
        "tumor", 
        "tumor en parafina (ffpe)",
        "tumor en parafina",
        "líquido pleural",
        "liquido pleural"
    ]:
        return "tumor"

    return "desconocida"

def clasificar_final_row(row):
    """
    Aplica la regla clínica combinada:
    Necesita:
      row["tipo_muestra_norm"]            -> sangre / tumor / desconocida
      row["es_reportada_somatica"]        -> True/False/NaN
    Devuelve (clasificacion_final, razon_clasificacion)
    """

    tm = row.get("tipo_muestra_norm", "desconocida")
    es_som = row.get("es_reportada_somatica", np.nan)

    # Caso 1: muestra tumor
    if tm == "tumor":
        if es_som is True:
            return "somática", "tumor + COSMIC somatic"
        else:
            # aunque COSMIC no la marque, al venir de tumor la mantenemos somática
            return "somática (solo por origen tumoral)", "tumor sin evidencia COSMIC somatic"

    # Caso 2: muestra sangre
    if tm == "sangre":
        if es_som is True:
            # sangre pero igual reportada somática en COSMIC
            return "somática", "sangre + COSMIC somatic"
        else:
            # sangre sin evidencia somática en COSMIC => germinal
            return "germinal", "sangre sin COSMIC somatic"

    # Caso 3: desconocida / otra
    if es_som is True:
        return "somática", "sin tipo_muestra claro + COSMIC somatic"

    return "indeterminado", "sin tipo_muestra claro y sin COSMIC somatic"


def main(path_variantes_extraidas, path_cosmic_tsv, path_out_csv):
    # ------------------------------------------------------------------
    # 1) Leer variantes_extraidas.csv
    # ------------------------------------------------------------------
    df = pd.read_csv(
        path_variantes_extraidas,
        na_values=["", "NA"],
        encoding="utf-8"
    )

    # chequeamos columnas clave
    required_cols = ["id", "variante"]
    for col in required_cols:
        if col not in df.columns:
            raise RuntimeError(f"El CSV de variantes no tiene la columna obligatoria '{col}'")

    # columna tipo_de_muestra (lo que tenés en tus datos reales)
    if "tipo_de_muestra" not in df.columns:
        raise RuntimeError(
            "El CSV de variantes no tiene columna 'tipo_de_muestra'. "
            "Ajustá el nombre acá si en tu base está distinto."
        )

    # ------------------------------------------------------------------
    # 2) Armar df_keys: parsear gene, hgvs_c, hgvs_p, aa3, aa1
    # ------------------------------------------------------------------
    df_keys = pd.DataFrame({
        "id":       df["id"],
        "variante": df["variante"],
        "gene":     df["variante"].apply(lambda x: extract_first(pat_gene, x)),
        "hgvs_c":   df["variante"].apply(lambda x: extract_first(pat_c, x)),
        "hgvs_p":   df["variante"].apply(lambda x: extract_first(pat_p, x)),
    })

    df_keys["aa3"] = df_keys["hgvs_p"].str.replace(r"^p\.", "", regex=True)
    df_keys["aa1"] = df_keys["aa3"].apply(aa3_to_aa1)

    # ------------------------------------------------------------------
    # 3) Subconjunto único con las llaves que consultaremos en COSMIC
    # ------------------------------------------------------------------
    keys_unique = (
        df_keys.loc[
            df_keys["gene"].notna() &
            (df_keys["aa3"].notna() | df_keys["aa1"].notna()),
            ["gene", "hgvs_c", "aa3", "aa1"]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    # ------------------------------------------------------------------
    # 4) DuckDB en memoria: cargar COSMIC TSV como vista
    # ------------------------------------------------------------------
    con = duckdb.connect(database=":memory:")

    con.execute("DROP VIEW IF EXISTS cosmic;")
    con.execute(f"""
        CREATE VIEW cosmic AS
        SELECT * FROM read_csv_auto('{path_cosmic_tsv}',
                                    delim='\t',
                                    header=true,
                                    sample_size=1000000);
    """)

    # subimos las keys para hacer join SQL
    con.register("keys_unique_df", keys_unique)
    con.execute("CREATE OR REPLACE TABLE t_input AS SELECT * FROM keys_unique_df;")

    # ------------------------------------------------------------------
    # 5) JOIN con COSMIC
    #    - detectamos si COSMIC marca la variante como somatic
    # ------------------------------------------------------------------
    qry = """
    WITH c_norm AS (
      SELECT
        GENE_SYMBOL, HGVSC, HGVSP, MUTATION_CDS, MUTATION_AA,
        MUTATION_ID, MUTATION_DESCRIPTION, MUTATION_SOMATIC_STATUS,
        REPLACE(split_part(HGVSP, ':', 2), 'p.', '') AS hgvsp_clean,
        REPLACE(MUTATION_AA, 'p.', '')               AS mutaa_clean
      FROM cosmic
    )
    SELECT
      t.gene, t.hgvs_c, t.aa3, t.aa1,
      COUNT(*) AS hits_total,
      SUM(CASE WHEN lower(c.MUTATION_SOMATIC_STATUS) LIKE '%somatic%' THEN 1 ELSE 0 END)
        AS hits_somatic
    FROM t_input t
    LEFT JOIN c_norm c
      ON  c.GENE_SYMBOL = t.gene
      AND (
            c.mutaa_clean = t.aa3
         OR c.mutaa_clean = t.aa1
         OR c.hgvsp_clean = t.aa3
         OR c.hgvsp_clean = t.aa1
         OR c.MUTATION_CDS = t.hgvs_c
         OR c.HGVSC LIKE ('%:' || t.hgvs_c)
         )
    GROUP BY 1,2,3,4
    ORDER BY hits_somatic DESC, hits_total DESC
    """

    cosmic_match = con.execute(qry).fetch_df()
    con.close()

    # ------------------------------------------------------------------
    # 6) Marcar "es_reportada_somatica" a nivel (gene,hgvs_c,aa3,aa1)
    # ------------------------------------------------------------------
    cosmic_flag = (
        cosmic_match
        .assign(es_reportada_somatica=lambda d: d["hits_somatic"] > 0)
        [["gene", "hgvs_c", "aa3", "aa1", "es_reportada_somatica"]]
    )

    # ------------------------------------------------------------------
    # 7) Volver a unir a df_keys (cada fila paciente-variante)
    # ------------------------------------------------------------------
    df_flag = df_keys.merge(
        cosmic_flag,
        on=["gene", "hgvs_c", "aa3", "aa1"],
        how="left"
    )

    # df_flag ahora tiene:
    # id, variante, gene, hgvs_c, aa3, aa1, es_reportada_somatica (True/False/NaN)

    # ------------------------------------------------------------------
    # 8) Reducir a una fila por id (porque un mismo id puede mapear
    #    a varias entradas parseadas).
    # ------------------------------------------------------------------
    df_flag_id = (
        df_flag.groupby("id", dropna=False)["es_reportada_somatica"]
        .apply(reduce_flag)  # <- True / False / NaN
        .reset_index()
    )

    # ------------------------------------------------------------------
    # 9) Unir eso de vuelta al dataframe completo original
    # ------------------------------------------------------------------
    df_final = df.merge(df_flag_id, on="id", how="left")

    # Normalizamos tipo_de_muestra → tipo_muestra_norm
    df_final["tipo_muestra_norm"] = df_final["tipo_de_muestra"].apply(normalizar_tipo_muestra)

    # ------------------------------------------------------------------
    # 10) Aplicar la regla clínica combinada fila a fila
    # ------------------------------------------------------------------
    clasifs = df_final.apply(clasificar_final_row, axis=1)
    df_final["clasificacion_final"] = clasifs.apply(lambda x: x[0])
    df_final["razon_clasificacion"] = clasifs.apply(lambda x: x[1])

    # ------------------------------------------------------------------
    # 11) Guardar resultado
    # ------------------------------------------------------------------
    df_final.to_csv(path_out_csv, index=False, encoding="utf-8")
    print(f"[OK] Archivo con clasificación COSMIC -> {path_out_csv}")
    print("Resumen clasificacion_final:")
    print(df_final["clasificacion_final"].value_counts(dropna=False))


if __name__ == "__main__":
    # Uso:
    #   python cosmic_classify.py <variantes_extraidas.csv> <cosmic.tsv> <salida.csv>
    if len(sys.argv) != 4:
        print("Uso: python cosmic_classify.py <variantes_extraidas.csv> <cosmic.tsv> <salida.csv>")
        sys.exit(1)

    variantes_csv = sys.argv[1]
    cosmic_tsv    = sys.argv[2]
    salida_csv    = sys.argv[3]

    main(variantes_csv, cosmic_tsv, salida_csv)
