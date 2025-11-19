#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
import math

def normalize_chrom_for_vcf(chrom_val):
    """
    Convierte:
      "17"     -> "chr17"
      "X"      -> "chrX"
      "Y"      -> "chrY"
      "MT"     -> "chrM"
      "M"      -> "chrM"
      "chr17"  -> "chr17" (lo dejamos igual)
    """
    if pd.isna(chrom_val):
        return None
    c = str(chrom_val).strip()

    # si ya empieza con "chr", lo dejamos
    if c.lower().startswith("chr"):
        base = c[3:]
    else:
        base = c

    # normalizamos mitocondrial
    if base in ["MT", "Mt", "mt", "M", "m", "MTDNA", "mtdna"]:
        base = "M"

    return f"chr{base}"

def pick_id(row):
    """
    Elegir ID de la variante para la columna ID del VCF.
    Preferimos Id_clinvar si existe/está no-vacío.
    Si no, fabricamos algo tipo var_<idInterno>
    """
    # Ajustá el nombre exacto si en tu CSV NO se llama 'Id_clinvar'
    clinvar_id = row.get("Id_clinvar", None)
    if pd.notna(clinvar_id) and str(clinvar_id).strip() != "":
        return str(clinvar_id)

    # fallback: usar id interno de la fila (observación)
    internal = row.get("id", None)
    if pd.notna(internal):
        return f"var_{internal}"

    # super fallback
    return "."

def build_info_field(group_df, source_label):
    """
    Arma el campo INFO para UNA variante ya agrupada.
    group_df = todas las filas (pacientes) que comparten misma coord REF/ALT.
    source_label = "SOMATIC" o "GERMINAL"

    INFO que metemos:
      SAMPLE_COUNT=...
      CLINVAR_MATCH=0/1 (1 si al menos una fila clinvar_match==True)
      CLINVAR_INFO=... (si coincidió ClinVar, tomamos el primer INFO_STRING no nulo)
      SOURCE=source_label
    """
    sample_count = len(group_df)

    clinvar_any = False
    clinvar_info_val = ""
    if "clinvar_match" in group_df.columns:
        # si alguna fila está matcheada
        clinvar_any = group_df["clinvar_match"].fillna(False).astype(bool).any()

    if clinvar_any and "INFO_STRING" in group_df.columns:
        # tomar la primera anotación ClinVar disponible
        first_info = group_df["INFO_STRING"].dropna()
        if len(first_info) > 0:
            clinvar_info_val = first_info.iloc[0]
        else:
            clinvar_info_val = ""
    else:
        clinvar_info_val = ""

    # sanitizar para que no rompa el VCF (sacar espacios/tab/nuevas líneas)
    # VCF INFO debe ser una sola línea sin tabs.
    clinvar_info_val = str(clinvar_info_val).replace("\t", " ").replace("\n", " ").replace("\r", " ")

    parts = [
        f"SAMPLE_COUNT={sample_count}",
        f"CLINVAR_MATCH={'1' if clinvar_any else '0'}",
        f"SOURCE={source_label}"
    ]
    if clinvar_info_val.strip() != "":
        # vamos a ponerlo como CLINVAR_INFO=...
        # IMPORTANTE: en VCF es buena práctica no meter '=' dentro del valor,
        # pero ClinVar mete 'ALLELEID=...' etc.
        # Lo dejamos crudo igual, porque esto es para screening interno.
        parts.append("CLINVAR_INFO=" + clinvar_info_val)

    return ";".join(parts)

def main(input_csv, source_label, output_vcf):
    """
    input_csv   = variantes_somaticas_API_anotadas.csv (o germinales)
    source_label= "SOMATIC" o "GERMINAL"
    output_vcf  = archivo .vcf a escribir (texto plano)
    """

    df = pd.read_csv(input_csv, dtype=str)  # leemos todo como string para evitar problemas de tipos
    # Convertir POS_from_hgvs a int para poder ordenar correctamente
    if "POS_from_hgvs" not in df.columns:
        raise RuntimeError("No encuentro la columna POS_from_hgvs en el CSV de entrada.")
    df["POS_from_hgvs"] = pd.to_numeric(df["POS_from_hgvs"], errors="coerce")

    # Normalizamos CHROM tipo "17" -> "chr17"
    if "CHROM_from_hgvs" not in df.columns:
        raise RuntimeError("No encuentro la columna CHROM_from_hgvs en el CSV de entrada.")
    df["CHROM_vcf"] = df["CHROM_from_hgvs"].apply(normalize_chrom_for_vcf)

    # Filtramos filas que tienen todo lo necesario
    needed_cols = ["CHROM_vcf","POS_from_hgvs","REF_from_hgvs","ALT_from_hgvs"]
    for c in needed_cols:
        if c not in df.columns:
            raise RuntimeError(f"Falta columna {c} en el CSV. No puedo construir VCF.")
    df_valid = df.dropna(subset=["CHROM_vcf","POS_from_hgvs","REF_from_hgvs","ALT_from_hgvs"])
    df_valid = df_valid[df_valid["REF_from_hgvs"].str.len() > 0]
    df_valid = df_valid[df_valid["ALT_from_hgvs"].str.len() > 0]

    # Agrupamos por variante única (CHROM,POS,REF,ALT)
    group_cols = ["CHROM_vcf","POS_from_hgvs","REF_from_hgvs","ALT_from_hgvs"]
    grouped = df_valid.groupby(group_cols, dropna=False)

    # Vamos a construir un DataFrame con columnas VCF:
    records = []
    for (chrom, pos, ref, alt), subdf in grouped:
        # ID preferente
        # usamos la primera fila del subdf como "representante" para ID
        rep = subdf.iloc[0]
        vcf_id = pick_id(rep)

        # QUAL
        vcf_qual = "."

        # FILTER
        vcf_filter = "PASS"

        # INFO
        vcf_info = build_info_field(subdf, source_label)

        # Guardar
        records.append({
            "CHROM":  chrom,
            "POS":    int(pos) if (not pd.isna(pos) and not isinstance(pos, str)) else int(float(pos)),
            "ID":     vcf_id,
            "REF":    str(ref),
            "ALT":    str(alt),
            "QUAL":   vcf_qual,
            "FILTER": vcf_filter,
            "INFO":   vcf_info
        })

    vcf_df = pd.DataFrame(records)

    # Orden VCF canónico: por cromosoma y posición
    # Intento de orden cromosómico humano: chr1..chr22,chrX,chrY,chrM
    def chrom_rank(c):
        c = c.replace("chr","")
        # manejamos M como 100
        if c == "X": return 23
        if c == "Y": return 24
        if c in ["M","MT"]: return 25
        # intentar número
        try:
            return int(c)
        except:
            return 99
    vcf_df["chrom_order"] = vcf_df["CHROM"].apply(chrom_rank)
    vcf_df = vcf_df.sort_values(by=["chrom_order","POS"]).drop(columns=["chrom_order"])

    # Escribimos VCF plano
    with open(output_vcf, "w", encoding="utf-8") as fh:
        # Header mínimo
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##reference=GRCh38\n")
        fh.write("##INFO=<ID=SAMPLE_COUNT,Number=1,Type=Integer,Description=\"How many patient-level rows carried this variant in this cohort\">\n")
        fh.write("##INFO=<ID=CLINVAR_MATCH,Number=1,Type=Integer,Description=\"1 if variant matched ClinVar entry in join_with_vcf step\">\n")
        fh.write("##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"SOMATIC or GERMINAL cohort label\">\n")
        fh.write("##INFO=<ID=CLINVAR_INFO,Number=1,Type=String,Description=\"ClinVar INFO_STRING subset for this variant (first hit)\">\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for _, row in vcf_df.iterrows():
            line = [
                row["CHROM"],
                str(row["POS"]),
                row["ID"],
                row["REF"],
                row["ALT"],
                row["QUAL"],
                row["FILTER"],
                row["INFO"]
            ]
            fh.write("\t".join(line) + "\n")

    print(f"[OK] VCF escrito en {output_vcf}")
    print(f"Total variantes únicas: {len(vcf_df)}")


if __name__ == "__main__":
    # Uso:
    #   python csv_to_vcf.py <input_csv> <SOMATIC|GERMINAL> <output.vcf>
    if len(sys.argv) != 4:
        print("Uso: python csv_to_vcf.py <input_csv> <SOMATIC|GERMINAL> <output.vcf>")
        sys.exit(1)

    in_csv = sys.argv[1]
    label  = sys.argv[2].upper()
    out_vcf = sys.argv[3]

    if label not in ["SOMATIC","GERMINAL"]:
        raise RuntimeError("El segundo argumento debe ser SOMATIC o GERMINAL")

    main(in_csv, label, out_vcf)
