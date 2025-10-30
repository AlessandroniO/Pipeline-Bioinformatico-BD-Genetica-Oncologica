import pandas as pd
import sys
import os

def analyze_missing_data(input_csv, output_report="data/results/qc/missing_data_report.csv"):
    """
    Analiza la completitud del archivo de variantes extraídas y genera un reporte.
    """

    print(f"📥 Leyendo archivo: {input_csv}")
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        print(f"❌ ERROR al leer el archivo CSV: {e}")
        sys.exit(1)

    print(f"✅ Archivo cargado con {len(df)} filas y {len(df.columns)} columnas.")

    # =========================================================
    # 1. Calcular métricas de completitud
    # =========================================================
    missing_counts = df.isnull().sum()
    total_rows = len(df)
    completeness = 100 * (1 - missing_counts / total_rows)

    report_df = pd.DataFrame({
        "Columna": df.columns,
        "Valores_nulos": missing_counts,
        "Porcentaje_completitud": completeness.round(2)
    }).sort_values(by="Porcentaje_completitud")

    # =========================================================
    # 2. Mostrar resumen en consola
    # =========================================================
    print("\n📊 Resumen de columnas con datos faltantes:")
    print(report_df.to_string(index=False))

    # =========================================================
    # 3. Guardar reporte como CSV
    # =========================================================
    os.makedirs(os.path.dirname(output_report), exist_ok=True)
    report_df.to_csv(output_report, index=False, encoding="utf-8-sig")
    print(f"\n💾 Reporte exportado a: {output_report}")

    # =========================================================
    # 4. Sugerencias automáticas
    # =========================================================
    print("\n🧩 Sugerencias iniciales:")
    for col in report_df["Columna"]:
        if col.lower() in ["chrom", "cromosoma"]:
            print(f" → Completar '{col}' desde ClinVar (campo CHROM o ID_variante).")
        elif col.lower() in ["pos", "posición"]:
            print(f" → Completar '{col}' desde ClinVar (campo POS).")
        elif col.lower() in ["ref", "alt"]:
            print(f" → Completar '{col}' con REF/ALT de ClinVar o simular bases aleatorias.")
        elif col.lower() == "tipo_de_muestra":
            print(f" → Usar tabla 'Muestra_genomica' o CIE-10 para clasificar como Germinal/Somática.")
        elif col.lower() in ["sexo", "edad"]:
            print(f" → Completar desde tabla 'Paciente' (campos Sexo y Fecha_nacimiento).")
        elif col.lower() == "barrio":
            print(f" → Completar desde tabla 'Direccion' asociada a cada paciente.")
        elif report_df.loc[report_df["Columna"] == col, "Porcentaje_completitud"].values[0] < 80:
            print(f" ⚠ '{col}' tiene menos del 80% de completitud — revisar fuente de datos.")

    print("\n✅ Análisis finalizado correctamente.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python analyze_missing_data.py <archivo_variantes_extraidas.csv>")
        sys.exit(1)

    input_csv = sys.argv[1]
    analyze_missing_data(input_csv)
