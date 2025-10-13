import sys
import pandas as pd
import os

# ===============================================================
# 1. Funciones auxiliares
# ===============================================================

def cargar_diccionario_cie10(cie10_path):
    """
    Carga la tabla CIE-10 local (Lista_tumores_CIE-10.csv)
    y genera un diccionario que mapea cada código a su tipo de mutación:
    somática o germinal.
    """
    cie10_df = pd.read_csv(cie10_path)
    cie10_df = cie10_df.rename(columns=str.lower)
    
    # Normalizamos los nombres esperados
    cols = cie10_df.columns
    codigo_col = [c for c in cols if "cie" in c or "codigo" in c][0]
    tipo_col = [c for c in cols if "tipo" in c][0]

    cie10_map = dict(zip(cie10_df[codigo_col], cie10_df[tipo_col].str.lower()))
    print(f"INFO: Diccionario CIE-10 cargado ({len(cie10_map)} códigos reconocidos).")
    return cie10_map


def inferir_tipo_mutacion(row, cie10_map):
    """
    Determina si una variante es somática o germinal según:
      1. Tipo_de_muestra
      2. Código CIE-10 (si disponible)
      3. Fuente (ClinVar, COSMIC, DECIPHER, etc.)
    """
    tipo_muestra = str(row.get("Tipo_de_muestra", "")).lower()
    cie10 = str(row.get("Codigo_cie10", "")).upper()
    fuente = str(row.get("Fuente", "")).lower()

    # --- Reglas jerárquicas ---
    # 1. Según fuente externa
    if "cosmic" in fuente:
        return "somatica"
    if "clinvar" in fuente or "decipher" in fuente or "spada" in fuente:
        return "germinal"

    # 2. Según tipo de muestra
    if "tumor" in tipo_muestra or "biopsia" in tipo_muestra:
        return "somatica"
    if "sangre" in tipo_muestra or "leucocito" in tipo_muestra:
        return "germinal"

    # 3. Según CIE-10
    if cie10 in cie10_map:
        return cie10_map[cie10]

    return "desconocida"


# ===============================================================
# 2. Lógica principal
# ===============================================================

def process_data():
    if len(sys.argv) < 4:
        print("Uso: python data_processing_split.py <input_csv> <output_dir> <cie10_file>")
        sys.exit(1)

    input_csv = sys.argv[1]
    output_dir = sys.argv[2]
    cie10_file = sys.argv[3]

    print(f"INFO: Leyendo datos de: {input_csv}")
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        print(f"❌ ERROR al leer el archivo CSV: {e}")
        sys.exit(1)

    initial_count = len(df)
    print(f"INFO: {initial_count} filas iniciales cargadas.")

    # --- Limpieza básica ---
    df = df.dropna(subset=["Id_muestra"], how="all")
    df["Tipo_de_muestra"] = df.get("Tipo_de_muestra", "").fillna("Desconocida")

    # --- Cargar diccionario CIE-10 ---
    cie10_map = cargar_diccionario_cie10(cie10_file)

    # --- Inferir tipo de mutación ---
    df["Tipo_mutacion"] = df.apply(lambda x: inferir_tipo_mutacion(x, cie10_map), axis=1)

    somatic_df = df[df["Tipo_mutacion"] == "somatica"]
    germline_df = df[df["Tipo_mutacion"] == "germinal"]
    unknown_df = df[df["Tipo_mutacion"] == "desconocida"]

    print(f"INFO: Somáticas: {len(somatic_df)} | Germinales: {len(germline_df)} | Desconocidas: {len(unknown_df)}")

    # --- Exportar resultados ---
    os.makedirs(output_dir, exist_ok=True)

    somatic_out = os.path.join(output_dir, "variantes_somaticas_limpias.csv")
    germline_out = os.path.join(output_dir, "variantes_germinales_limpias.csv")
    unknown_out = os.path.join(output_dir, "variantes_tipo_desconocido.csv")

    somatic_df.to_csv(somatic_out, index=False)
    germline_df.to_csv(germline_out, index=False)
    unknown_df.to_csv(unknown_out, index=False)

    print(f"✅ Exportados:\n  - {somatic_out}\n  - {germline_out}\n  - {unknown_out}")

# ===============================================================
# 3. Ejecución directa
# ===============================================================

if __name__ == "__main__":
    process_data()
