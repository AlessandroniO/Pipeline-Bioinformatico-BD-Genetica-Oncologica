# cosmic_classify.py
# Requiere: pandas, duckdb, pyarrow, numpy, regex
# Objetivo:
#   - Leer variantes_extraidas.csv
#   - Cruzar con COSMIC (archivo TSV grande)
#   - Agregar columnas:
#       es_reportada_somatica (True/False/NaN)
#       tipo_celular_afectado ("somática" / "no somática reportada" / "indeterminado")
#   - Guardar CSV final enriquecido

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
    series es la columna es_reportada_somatica para un mismo id.
    - Si hay al menos un True -> True
    - Si TODO es NaN          -> NaN
    - Sino                    -> False
    """
    if series.any(skipna=True):
        return True
    if series.isna().all():
        return np.nan
    return False


def main(path_variantes_extraidas, path_cosmic_tsv, path_out_csv):
    # 1) Leer variantes_extraidas.csv
    df = pd.read_csv(
        path_variantes_extraidas,
        na_values=["", "NA"],
        encoding="utf-8"
    )

    # 2) Armar df_keys con gene, hgvs_c, hgvs_p, aa3, aa1
    df_keys = pd.DataFrame({
        "id":       df["id"],
        "variante": df["variante"],
        "gene":     df["variante"].apply(lambda x: extract_first(pat_gene, x)),
        "hgvs_c":   df["variante"].apply(lambda x: extract_first(pat_c, x)),
        "hgvs_p":   df["variante"].apply(lambda x: extract_first(pat_p, x)),
    })

    df_keys["aa3"] = df_keys["hgvs_p"].str.replace(r"^p\.", "", regex=True)
    df_keys["aa1"] = df_keys["aa3"].apply(aa3_to_aa1)

    # 3) Subconjunto único de llaves que tengan gene y (aa3 o aa1)
    keys_unique = (
        df_keys.loc[
            df_keys["gene"].notna() &
            (df_keys["aa3"].notna() | df_keys["aa1"].notna()),
            ["gene", "hgvs_c", "aa3", "aa1"]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    # 4) DuckDB: cargar COSMIC TSV como vista
    # IMPORTANTE: no creamos base persistente, usamos in-memory
    con = duckdb.connect(database=":memory:")

    # Creamos VIEW 'cosmic' leyendo el TSV que vos ya descargaste
    con.execute("DROP VIEW IF EXISTS cosmic;")
    con.execute(f"""
        CREATE VIEW cosmic AS
        SELECT * FROM read_csv_auto('{path_cosmic_tsv}', delim='\t', header=true, sample_size=1000000);
    """)

    # Registramos las keys únicas como tabla temporal t_input
    con.register("keys_unique_df", keys_unique)
    con.execute("CREATE OR REPLACE TABLE t_input AS SELECT * FROM keys_unique_df;")

    # 5) JOIN lógica tipo R
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

    # 6) Marcar si aparece como somática en COSMIC
    cosmic_flag = (
        cosmic_match
        .assign(es_reportada_somatica=lambda d: d["hits_somatic"] > 0)
        [["gene", "hgvs_c", "aa3", "aa1", "es_reportada_somatica"]]
    )

    # 7) Volver a unir a df_keys fila a fila
    df_flag = df_keys.merge(
        cosmic_flag,
        on=["gene", "hgvs_c", "aa3", "aa1"],
        how="left"
    )

    # 8) Reducir a "por id" (id = tu identificador de variante en esa tabla)
    df_flag_id = (
        df_flag.groupby("id", dropna=False)["es_reportada_somatica"]
        .apply(reduce_flag)
        .reset_index()
    )

    # 9) Volcar esa info en el data frame original
    df_final = df.merge(df_flag_id, on="id", how="left")

    df_final["tipo_celular_afectado"] = np.select(
        [
            df_final["es_reportada_somatica"] == True,
            df_final["es_reportada_somatica"] == False,
            df_final["es_reportada_somatica"].isna()
        ],
        ["somática", "no somática reportada", "indeterminado"],
        default="indeterminado"
    )

    # 10) Guardar resultado
    df_final.to_csv(path_out_csv, index=False, encoding="utf-8")

    con.close()
    print(f"[OK] Archivo con clasificación COSMIC -> {path_out_csv}")


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
