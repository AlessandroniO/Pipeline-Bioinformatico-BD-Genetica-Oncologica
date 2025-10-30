import sys
import pandas as pd
import re


def explode_hgvs_g_to_coords(hgvs_g: str):
    """
    Ejemplo entrada:
        'NC_000018.10:g.63637928C>T'
    Salida:
        chrom='18', pos=63637928, ref='C', alt='T'
    Para variantes no-SNV simples devolvemos (None,None,None,None).
    """
    if not isinstance(hgvs_g, str):
        return None, None, None, None

    if not (hgvs_g.startswith("NC_") and ":g." in hgvs_g and ">" in hgvs_g):
        return None, None, None, None

    left, change = hgvs_g.split(":g.", 1)

    m = re.match(r"(\d+)([ACGT])>([ACGT])", change)
    if not m:
        return None, None, None, None

    pos = int(m.group(1))
    ref = m.group(2)
    alt = m.group(3)

    chrom_raw = left  # ej "NC_000018.10"
    chrom_num = chrom_raw.replace("NC_", "").split(".")[0]  # "000018"
    chrom_num = chrom_num.lstrip("0")                       # "18"
    if chrom_num == "":
        chrom_num = "0"

    return chrom_num, pos, ref, alt


def normalize_chrom(x: str):
    """
    Fuerza cromosoma a formato comparable:
    - '18' -> '18'
    - 'chr18' -> '18'
    - 'NC_000018.10' -> '18'
    - '000018' -> '18'
    - 'chrX' -> 'X'
    """
    if x is None:
        return None
    s = str(x).strip()
    if s == "":
        return None

    # quitar 'chr'
    if s.lower().startswith("chr"):
        s = s[3:]

    # NC_000018.10 -> 000018
    if s.startswith("NC_"):
        core = s.replace("NC_", "").split(".")[0]
        s = core

    # sacar ceros a la izquierda si es numérico
    if s.isdigit():
        s_no0 = s.lstrip("0")
        s = s_no0 if s_no0 != "" else "0"

    # cromosoma mitocondrial, variantes posibles
    if s.upper() in ["M", "MT", "MTDNA"]:
        s = "MT"

    return s


def safe_pos(series):
    """
    Convierte posiciones a enteros pandas nullable (Int64).
    """
    return pd.to_numeric(series, errors="coerce").astype("Int64")


def main(sample_csv, clinvar_tsv, output_csv):
    # Leemos tus variantes (fila = hallazgo en un paciente)
    df_sample = pd.read_csv(sample_csv)

    # Leemos ClinVar subset
    df_clinvar = pd.read_csv(clinvar_tsv, sep="\t", dtype=str)

    required_cols = ["#CHROM", "POS", "REF", "ALT", "INFO_STRING"]
    for c in required_cols:
        if c not in df_clinvar.columns:
            raise RuntimeError(f"En ClinVar TSV falta la columna {c}")

    # Asegurar tipos básicos en ClinVar
    df_clinvar["POS"] = safe_pos(df_clinvar["POS"])
    df_clinvar["REF"] = df_clinvar["REF"].astype(str)
    df_clinvar["ALT"] = df_clinvar["ALT"].astype(str)

    # ----------------------------------------------------
    # 1. Derivar coordenadas genómicas (chrom,pos,ref,alt)
    #    para CADA FILA de df_sample si aún no existen.
    # ----------------------------------------------------
    if not {"CHROM_from_hgvs","POS_from_hgvs","REF_from_hgvs","ALT_from_hgvs"}.issubset(df_sample.columns):
        chroms, poss, refs, alts = [], [], [], []
        for g in df_sample["HGVS_G_GRCh38"]:
            c, p, r, a = explode_hgvs_g_to_coords(g)
            chroms.append(c)
            poss.append(p)
            refs.append(r)
            alts.append(a)
        df_sample["CHROM_from_hgvs"] = chroms
        df_sample["POS_from_hgvs"]   = poss
        df_sample["REF_from_hgvs"]   = refs
        df_sample["ALT_from_hgvs"]   = alts

    # Normalizamos cromosoma y posición en df_sample
    df_sample["CHROM_norm"] = df_sample["CHROM_from_hgvs"].apply(normalize_chrom)
    df_sample["POS_norm"]   = safe_pos(df_sample["POS_from_hgvs"])
    df_sample["REF_from_hgvs"] = df_sample["REF_from_hgvs"].astype(str)
    df_sample["ALT_from_hgvs"] = df_sample["ALT_from_hgvs"].astype(str)

    # ----------------------------------------------------
    # 2. Creamos tabla de variantes ÚNICAS por coordenada
    #    clave = (CHROM_norm, POS_norm, REF_from_hgvs, ALT_from_hgvs)
    # ----------------------------------------------------
    unique_vars = (
        df_sample[
            ["CHROM_norm","POS_norm","REF_from_hgvs","ALT_from_hgvs"]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    # ----------------------------------------------------
    # 3. Preparamos ClinVar con las MISMAS llaves normalizadas
    #    - Generamos CHROM_norm en ClinVar
    # ----------------------------------------------------
    df_clinvar["CHROM_norm"] = df_clinvar["#CHROM"].apply(normalize_chrom)

    # Ya tenemos:
    #   df_clinvar["POS"]  (Int64)
    #   df_clinvar["REF"]  (str)
    #   df_clinvar["ALT"]  (str)

    # ----------------------------------------------------
    # 4. Hacemos el merge SOLO UNA VEZ:
    #    unique_vars  ⨝  df_clinvar
    # ----------------------------------------------------
    unique_annot = pd.merge(
        unique_vars,
        df_clinvar,
        how="left",
        left_on=["CHROM_norm","POS_norm","REF_from_hgvs","ALT_from_hgvs"],
        right_on=["CHROM_norm","POS","REF","ALT"]
    )

    # Agregamos flag clinvar_match en esta tabla reducida
    unique_annot["clinvar_match"] = unique_annot["INFO_STRING"].notna()

    # ----------------------------------------------------
    # 5. Propagamos la anotación de unique_annot
    #    de vuelta a CADA FILA de df_sample.
    #
    #    Hacemos un merge otra vez, pero ahora entre:
    #        df_sample (todas las filas/pacientes)
    #    y   unique_annot (cada variante única con INFO_STRING)
    #
    #    Este merge es mucho más chico que hacerlo directo
    #    contra todo ClinVar fila a fila.
    # ----------------------------------------------------
    merged_full = pd.merge(
        df_sample,
        unique_annot[
            [
                "CHROM_norm",
                "POS_norm",
                "REF_from_hgvs",
                "ALT_from_hgvs",
                "INFO_STRING",
                "clinvar_match",
                "#CHROM",
                "POS",
                "REF",
                "ALT"
            ]
        ],
        how="left",
        on=["CHROM_norm","POS_norm","REF_from_hgvs","ALT_from_hgvs"]
    )

    # ----------------------------------------------------
    # 6. Guardar resultado final
    # ----------------------------------------------------
    merged_full.to_csv(output_csv, index=False)
    print(f"Anotación ClinVar lista -> {output_csv}")
    print(f"Total filas originales       : {len(df_sample)}")
    print(f"Variantes únicas anotadas    : {len(unique_vars)}")
    print(f"De esas, mapeadas en ClinVar : {unique_annot['clinvar_match'].sum()} (clinvar_match=True)")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python join_with_vcf.py <variantes_*_genomic.csv> <clinvar_grch38_parsed.tsv> <salida.csv>")
        sys.exit(1)

    sample_csv = sys.argv[1]
    clinvar_tsv = sys.argv[2]
    out_csv     = sys.argv[3]
    main(sample_csv, clinvar_tsv, out_csv)
