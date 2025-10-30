#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
import os
import csv

# ----------------------------
# Aux: Leer CSV con detección de separador y pruebas de encoding
# ----------------------------
def robust_read_csv(path, encodings=("utf-8","latin-1","cp1252")):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Archivo no encontrado: {path}")
    for enc in encodings:
        try:
            with open(path, "r", encoding=enc, errors="replace") as fh:
                sample = fh.read(8192)
                try:
                    dialect = csv.Sniffer().sniff(sample, delimiters=[",",";","\t","|"])
                    sep = dialect.delimiter
                except Exception:
                    sep = ","
            df = pd.read_csv(path, encoding=enc, sep=sep, low_memory=False)
            print(f"INFO: Leído {path} con encoding={enc} sep='{sep}' shape={df.shape}")
            return df
        except Exception as e:
            last_exc = e
    raise last_exc

# ----------------------------
# CIE-10 loader
# ----------------------------
def cargar_diccionario_cie10(cie10_path):
    if not os.path.exists(cie10_path):
        print(f"WARNING: Archivo CIE-10 no encontrado: {cie10_path}")
        return {}
    try:
        cie10_df = robust_read_csv(cie10_path)
        cie10_df.columns = cie10_df.columns.str.lower()
        codigo_col = next((c for c in cie10_df.columns if "cie" in c or "codigo" in c), None)
        tipo_col = next((c for c in cie10_df.columns if "tipo" in c), None)
        if not codigo_col or not tipo_col:
            print("WARNING: Columnas esperadas en CIE-10 no encontradas.")
            return {}
        cie10_map = dict(zip(cie10_df[codigo_col], cie10_df[tipo_col].astype(str).str.lower()))
        print(f"INFO: Diccionario CIE-10 cargado ({len(cie10_map)} códigos).")
        return cie10_map
    except Exception as e:
        print(f"ERROR al cargar CIE-10: {e}")
        return {}

# ----------------------------
# Inferencia tipo mutacion
# ----------------------------
def inferir_tipo_mutacion(row, cie10_map):
    tipo_muestra = str(row.get("tipo_de_muestra", "")).strip().lower()
    tipo_tumor = str(row.get("tipo_tumor", "")).strip().lower()
    mapeo_muestra = {
        "sangre": "germinal",
        "plasma": "germinal",
        "suero": "germinal",
        "liquido pleural": "somatica",
        "tumor en parafina": "somatica",
        "tumor": "somatica",
    }
    if tipo_muestra in mapeo_muestra:
        return mapeo_muestra[tipo_muestra]
    if tipo_tumor and tipo_tumor not in ["nan","desconocido",""]:
        return "somatica"
    return "desconocida"

# ----------------------------
# Main processing
# ----------------------------
def process_data():
    if len(sys.argv) != 6:
        print("Uso: python process_data.py <db_csv> <vcf_csv> <cie10_file> <somatic_out> <germline_out>")
        sys.exit(1)

    db_path, vcf_path, cie10_path, somatic_out, germline_out = sys.argv[1:]

    # Cargar CSVs robustamente
    try:
        df_db = robust_read_csv(db_path)
        df_vcf = robust_read_csv(vcf_path)
    except Exception as e:
        print(f"❌ ERROR al leer CSVs: {e}")
        sys.exit(1)

    df_db.columns = df_db.columns.str.lower()
    df_vcf.columns = df_vcf.columns.str.lower()
    print(f"INFO: Filas BD={len(df_db)} | Filas VCF={len(df_vcf)}")

    # QC columnas y nulos
    os.makedirs("data/results/qc", exist_ok=True)
    pd.DataFrame({"columns_db": df_db.columns}).to_csv("data/results/qc/db_columns.csv", index=False)
    pd.DataFrame({"columns_vcf": df_vcf.columns}).to_csv("data/results/qc/vcf_columns.csv", index=False)
    df_db.isna().sum().to_frame("n_missing").to_csv("data/results/qc/db_missing_counts.csv")
    df_vcf.isna().sum().to_frame("n_missing").to_csv("data/results/qc/vcf_missing_counts.csv")

    # --------------------- MERGE ----------------------
    merge_key = ['chrom', 'pos', 'ref', 'alt']
    for key in merge_key:
        if key not in df_db.columns:
            df_db[key] = pd.NA
        if key not in df_vcf.columns:
            df_vcf[key] = pd.NA

    # Renombrar columnas VCF para evitar colisiones
    df_vcf = df_vcf.rename(columns={
        'qual':'vcf_qual','filter':'vcf_filter','id':'vcf_id','tipo_variante':'vcf_tipo_variante'
    })

    # Merge left: mantener BD y agregar info ClinVar
    df_merged = pd.merge(df_db, df_vcf, on=merge_key, how='left')

    # Rellenar campos críticos desde VCF si faltan
    fill_map = {
        'qual':'vcf_qual',
        'filter':'vcf_filter',
        'id':'vcf_id',
        'tipo_variante':'vcf_tipo_variante'
    }
    for target, source in fill_map.items():
        if target in df_merged.columns and source in df_merged.columns:
            df_merged[target] = df_merged[target].fillna(df_merged[source])

    # Eliminar columnas auxiliares
    drop_cols = ['vcf_qual','vcf_filter','vcf_id','vcf_tipo_variante']
    df_merged.drop(columns=[c for c in drop_cols if c in df_merged.columns], inplace=True)

    # --------------------- CLASIFICACIÓN ----------------------
    cie10_map = cargar_diccionario_cie10(cie10_path)
    df_merged["tipo_mutacion"] = df_merged.apply(lambda x: inferir_tipo_mutacion(x, cie10_map), axis=1)

    # --------------------- CONTROL DE CALIDAD ----------------------
    qc_path = "data/results/qc/extraction_nulls.csv"
    df_merged["flag_null_qc"] = False
    if os.path.exists(qc_path):
        try:
            qc_df = robust_read_csv(qc_path)
            qc_df.columns = qc_df.columns.str.lower()
            qc_rows = df_merged.merge(qc_df, on=merge_key, how='inner')
            df_merged.loc[df_merged.index.isin(qc_rows.index), "flag_null_qc"] = True
            print(f"QC: {qc_rows.shape[0]} filas marcadas flag_null_qc=True")
        except Exception as e:
            print(f"ERROR QC integration: {e}")

    # --------------------- DIVISIÓN ----------------------
    somatic_df = df_merged[df_merged["tipo_mutacion"]=="somatica"]
    germline_df = df_merged[df_merged["tipo_mutacion"]=="germinal"]

    os.makedirs(os.path.dirname(somatic_out) or ".", exist_ok=True)
    somatic_df.to_csv(somatic_out, index=False)
    germline_df.to_csv(germline_out, index=False)

    # --------------------- RESUMEN QC ----------------------
    summary_list = []
    for tipo in ["somatica","germinal","desconocida"]:
        subset = df_merged[df_merged["tipo_mutacion"]==tipo]
        if subset.empty:
            continue
        total = len(subset)
        flagged = int(subset["flag_null_qc"].sum())
        perc_flagged = round(flagged/total*100,2) if total>0 else 0.0
        qc_status = "Revisión" if perc_flagged>5 else "OK"
        mean_qual = round(subset["qual"].mean(skipna=True),2) if "qual" in subset else None
        mean_frec = round(subset["frec_alelica"].mean(skipna=True),4) if "frec_alelica" in subset else None
        summary_list.append({
            "Tipo_mutacion":tipo,
            "Flagged":flagged,
            "Total":total,
            "Porcentaje_Flagged":perc_flagged,
            "Media_QUAL":mean_qual,
            "Frecuencia_Alelica_Promedio":mean_frec,
            "QC_Status":qc_status
        })

    qc_summary = pd.DataFrame(summary_list)
    qc_summary_path = "data/results/qc/qc_summary.csv"
    os.makedirs(os.path.dirname(qc_summary_path) or ".", exist_ok=True)
    qc_summary.to_csv(qc_summary_path, index=False)

    print(f"✅ Exportados: Somáticas={len(somatic_df)} -> {somatic_out}, Germinales={len(germline_df)} -> {germline_out}")
    print(f"📊 Resumen QC: {qc_summary_path}")

if __name__ == "__main__":
    process_data()
