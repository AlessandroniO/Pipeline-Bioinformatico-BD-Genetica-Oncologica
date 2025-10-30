#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import pandas as pd
import os
import re # Necesario para la corrección del CIE-10

# ===============================================================
# 1. Funciones auxiliares
# ===============================================================

def cargar_diccionario_cie10(cie10_path):
    """Carga la tabla CIE-10 y genera un diccionario {codigo: Ddescripcion}."""
    # NOTA: Esta función se mantiene para evitar fallos, aunque no se usa en inferir_tipo_mutacion.
    if not os.path.exists(cie10_path):
        print(f"❌ ERROR CRÍTICO: Archivo CIE-10 no encontrado en {cie10_path}")
        return {}

    try:
        # Intentar leer con UTF-8 por si el archivo fue guardado con esa codificación
        cie10_df = pd.read_csv(cie10_path, encoding='utf-8') 
        cie10_df.columns = cie10_df.columns.str.lower()
        
        # CORRECCIÓN 1: El error NameError se ocultaba tras el error list index out of range.
        # Se asegura que la columna de código existe.
        codigo_cols = [c for c in cie10_df.columns if "cie" in c or "codigo" in c]
        if not codigo_cols:
             raise ValueError("No se encontró columna 'Codigo'.")
        codigo_col = codigo_cols[0]

        # CORRECCIÓN 2: Uso de expresión regular para encontrar 'Ddescripcion'
        # El patrón 'Ddesc' que usaste busca la D mayúscula. Hay que buscar minúsculas.
        desc_cols = [c for c in cie10_df.columns if re.search(r'desc', c)]
        if not desc_cols:
             raise ValueError("No se encontró columna con 'desc' (ej. 'ddescripcion').")
        tipo_col = desc_cols[0] 

        cie10_map = dict(zip(cie10_df[codigo_col], cie10_df[tipo_col].astype(str).str.lower()))
        print(f"INFO: Diccionario CIE-10 cargado ({len(cie10_map)} códigos reconocidos).")
        return cie10_map
    except Exception as e:
        print(f"❌ ERROR al cargar o procesar el archivo CIE-10: {e}")
        # Se devuelve un diccionario vacío para que el script pueda continuar.
        return {} 

def inferir_tipo_mutacion(row, cie10_map):
    """Determina si una variante es somática o germinal según tipo_de_muestra."""
    
    # 1. DEFINICIÓN DE VARIABLE (CORRECCIÓN CRÍTICA DE NameError)
    # Se define la variable tipo_muestra DENTRO de esta función,
    # usando .get() para ser seguro, y se normaliza.
    tipo_muestra = str(row.get("tipo_de_muestra", "")).strip().lower()

    # 2. Mapeo de muestras (Incluye correcciones de formato y codificación)
    mapeo_muestra = {
        "sangre": "germinal",
        "liquido pleural": "somatica",
        "lÃ­quido pleural": "somatica",              # Corrección de Codificación (Ã­)
        "tumor en parafina (ffpe)": "somatica",     # Corregido el formato
        "tumor en parafina": "somatica",
        "tumor": "somatica",
    }
    
    # --------------------- CLASIFICACIÓN ---------------------
    
    # Prioridad 1: Tipo de Muestra
    if tipo_muestra in mapeo_muestra:
        return mapeo_muestra[tipo_muestra]
    
    # 3. Se elimina completamente la lógica de tipo_tumor.
    
    return "desconocida"

# ===============================================================
# 2. Lógica principal (SOLO PROCESAMIENTO)
# ===============================================================

def process_data():
    # Verificación de argumentos
    if len(sys.argv) != 5: 
        print("Uso: python data_processing.py <db_csv> <cie10_file> <somatic_out> <germline_out>")
        sys.exit(1)

    db_data_path, cie10_path, somatic_out, germline_out = sys.argv[1:]

    # Lectura de CSV de la base de datos
    try:
        # Se usa solo un pd.read_csv y se intenta forzar la codificación
        df_merged = pd.read_csv(db_data_path, encoding='utf-8') 
    except Exception as e:
        print(f"❌ ERROR al leer archivo DB CSV: {e}")
        sys.exit(1)

    df_merged.columns = df_merged.columns.str.lower()

    print(f"INFO: Filas de la Base de Datos cargadas: {len(df_merged)}")
    
    # --------------------- CLASIFICACIÓN ----------------------
    cie10_map = cargar_diccionario_cie10(cie10_path) # Se carga el diccionario

    # La función .apply() utiliza la función auxiliar para llenar la columna 'tipo_mutacion'.
    df_merged["tipo_mutacion"] = df_merged.apply(lambda row: inferir_tipo_mutacion(row, cie10_map), axis=1)

    # --------------------- CONTROL DE CALIDAD (SIMPLIFICADO) ----------------------
    df_merged["flag_null_qc"] = False
    print("INFO: Se omitió la integración de control de calidad externo. 'flag_null_qc' es False por defecto.")

    # --------------------- DIVISIÓN ----------------------
    somatic_df = df_merged[df_merged["tipo_mutacion"] == "somatica"]
    germline_df = df_merged[df_merged["tipo_mutacion"] == "germinal"]
    unknown_df = df_merged[df_merged["tipo_mutacion"] == "desconocida"]

    os.makedirs(os.path.dirname(somatic_out), exist_ok=True)
    somatic_df.to_csv(somatic_out, index=False)
    germline_df.to_csv(germline_out, index=False)

    # --------------------- RESUMEN QC EXTENDIDO ----------------------
    # ... (El resto del resumen QC se mantiene, operando sobre somatic_df, germline_df, etc.) ...
    
    summary_list = []
    for tipo in ["somatica", "germinal", "desconocida"]:
        subset = df_merged[df_merged["tipo_mutacion"] == tipo]
        if subset.empty:
            continue
        total = len(subset)
        
        flagged = 0
        no_flagged = total 
        perc_flagged = 0.00
        qc_status = "OK" 

        # Aseguramos que los campos existan antes de calcular la media
        mean_qual = round(subset["qual"].mean(skipna=True), 2) if "qual" in subset else None
        mean_frec = round(subset["frec_alelica"].mean(skipna=True), 4) if "frec_alelica" in subset else None

        summary_list.append({
            "Tipo_mutacion": tipo,
            "No_Flagged": no_flagged,
            "Flagged": flagged,
            "Total": total,
            "Porcentaje_Flagged": perc_flagged,
            "Media_QUAL": mean_qual,
            "Frecuencia_Alelica_Promedio": mean_frec,
            "QC_Status": qc_status
        })

    qc_summary = pd.DataFrame(summary_list)
    qc_summary_path = "data/results/qc/qc_summary.csv"
    os.makedirs(os.path.dirname(qc_summary_path), exist_ok=True)
    qc_summary.to_csv(qc_summary_path, index=False)

    print(f"""
 Exportados:
  - Somáticas: {len(somatic_df)} → {somatic_out}
  - Germinales: {len(germline_df)} → {germline_out}
  - Desconocidas: {len(unknown_df)}
  - Filas marcadas con flag_null_qc: {df_merged['flag_null_qc'].sum()}

Resumen QC extendido guardado en:
  → {qc_summary_path}
""")

# ===============================================================
# 3. Ejecución directa
# ===============================================================

if __name__ == "__main__":
    process_data()