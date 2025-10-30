#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import pandas as pd

# =========================================================================
# 1. FUNCIÓN PRINCIPAL
# =========================================================================
def create_vcf(input_csv_path, output_vcf_path):
    """
    Convierte un CSV de variantes (limpio y QC pasado) a VCF estándar v4.3.
    """
    # Crear directorio de salida
    output_dir = os.path.dirname(output_vcf_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Leer CSV
    try:
        df = pd.read_csv(input_csv_path)
        print(f"INFO: Leyendo datos limpios de {input_csv_path} ({len(df)} filas)")
    except Exception as e:
        print(f"❌ ERROR al leer CSV: {e}")
        sys.exit(1)

    # Verificación de columnas críticas
    required_cols = ['chrom','pos','id_variante','ref','alt','tipo_tumor','frec_alelica','id_muestra']
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        print(f"❌ ERROR: Faltan columnas críticas en CSV: {missing_cols}")
        sys.exit(1)

    try:
        with open(output_vcf_path, 'w') as vcf_file:
            # --- Encabezado VCF ---
            vcf_file.write("##fileformat=VCFv4.3\n")
            vcf_file.write("##source=DW_Extraction_Pipeline\n")
            vcf_file.write('##INFO=<ID=DBID,Number=1,Type=String,Description="ID de la variante en la BD">\n')
            vcf_file.write('##INFO=<ID=TUMOR,Number=1,Type=String,Description="Tipo de tumor asociado">\n')
            vcf_file.write('##INFO=<ID=CS,Number=1,Type=String,Description="Significancia clínica">\n')
            vcf_file.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Frecuencia alélica">\n')
            vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

            # --- Cuerpo VCF ---
            for _, row in df.iterrows():
                chrom = str(row['chrom']).replace('chr','')
                pos = int(row['pos'])
                ref = str(row['ref']) if pd.notna(row['ref']) else "N"
                alt = str(row['alt']) if pd.notna(row['alt']) else "V"
                vcf_id = str(row['id_variante'])
                tumor = str(row['tipo_tumor']) if pd.notna(row['tipo_tumor']) else "ND"
                af_val = float(row['frec_alelica']) if pd.notna(row['frec_alelica']) else 0.0
                sample = str(row['id_muestra'])

                info_list = [
                    f"DBID={vcf_id}",
                    f"TUMOR={tumor}",
                    f"CS=ND",
                    f"AF={af_val:.4f}"
                ]
                info_str = ";".join(info_list)
                format_str = "GT:AF"
                sample_str = f"1/1:{af_val:.4f}"

                vcf_line = [
                    chrom,
                    pos,
                    vcf_id,
                    ref,
                    alt,
                    ".",          # QUAL
                    "PASS",       # FILTER
                    info_str,
                    format_str,
                    sample_str
                ]

                vcf_file.write("\t".join(map(str, vcf_line)) + "\n")

        print(f"✅ Éxito: VCF generado en {output_vcf_path} ({len(df)} variantes)")

    except Exception as e:
        print(f"❌ ERROR durante la conversión a VCF: {e}")
        sys.exit(1)

# =========================================================================
# 2. EJECUCIÓN
# =========================================================================
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python vcf_conversion.py <input_csv> <output_vcf>")
        sys.exit(1)
    input_csv = sys.argv[1]
    output_vcf = sys.argv[2]
    create_vcf(input_csv, output_vcf)
