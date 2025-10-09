import sys
import os
import pandas as pd
import numpy as np

# =========================================================================
# 1. FUNCIÓN PRINCIPAL DE PROCESAMIENTO
# =========================================================================
def process_and_split_data(input_csv_path, output_somatic_path, output_germline_path):
    """
    Lee el CSV extraído de la base de datos, limpia las cabeceras, 
    normaliza el tipo de muestra y divide el DataFrame en somático y germinal.
    
    Args:
        input_csv_path (str): Ruta al CSV extraído.
        output_somatic_path (str): Ruta para guardar el CSV de variantes somáticas.
        output_germline_path (str): Ruta para guardar el CSV de variantes germinales.
    """
    
    # 1. Crear directorios de salida
    os.makedirs(os.path.dirname(output_somatic_path), exist_ok=True)

    try:
        # 2. Leer los datos
        print(f"INFO: Leyendo datos desde {input_csv_path}...")
        df = pd.read_csv(input_csv_path)
        print(f"INFO: Total de {len(df)} filas leídas.")

    except FileNotFoundError:
        print(f"❌ ERROR CRÍTICO: Archivo de entrada no encontrado en {input_csv_path}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ ERROR CRÍTICO al leer el CSV: {e}")
        sys.exit(1)

    # 3. CORRECCIÓN CRUCIAL: Normalizar todas las cabeceras a minúsculas
    df.columns = df.columns.str.lower()
    
    # 4. Verificar la columna clave de división
    key_col = 'tipo_de_muestra'
    if key_col not in df.columns:
        print(f"❌ ERROR CRÍTICO: La columna '{key_col}' no fue encontrada en los datos.")
        print(f"Columnas encontradas: {list(df.columns)}")
        print("Esto indica un fallo en la extracción (db_extraction.py).")
        sys.exit(1)
        
    # 5. Normalizar los valores de la columna clave
    # Esto asegura que "SOMATICO", "Somático" o "somático " se traten como "somatico"
    df['muestra_lower'] = df[key_col].astype(str).str.strip().str.lower()
    
    # 6. Definir las categorías de división
    # La extracción debe garantizar que solo existan estos tipos o se debe 
    # añadir más lógica de limpieza.
    SOMATIC_KEYWORDS = ['somatica', 'somatico', 'tumor']
    GERMLINE_KEYWORDS = ['germinal', 'normal', 'sangre']
    
    # 7. Dividir los datos
    
    # Variables de salida
    somatic_data = df[df['muestra_lower'].isin(SOMATIC_KEYWORDS)]
    germline_data = df[df['muestra_lower'].isin(GERMLINE_KEYWORDS)]
    
    # 8. Guardar los archivos de salida
    
    # Guardar somáticas
    somatic_data.drop(columns=['muestra_lower']).to_csv(output_somatic_path, index=False)
    print(f"✅ Éxito: {len(somatic_data)} variantes somáticas guardadas en {output_somatic_path}")
    
    # Guardar germinales
    germline_data.drop(columns=['muestra_lower']).to_csv(output_germline_path, index=False)
    print(f"✅ Éxito: {len(germline_data)} variantes germinales guardadas en {output_germline_path}")

# =========================================================================
# 2. EJECUCIÓN
# =========================================================================
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python data_processing.py <input_csv> <output_somatic_csv> <output_germline_csv>")
        sys.exit(1)
    
    input_csv = sys.argv[1]
    output_somatic = sys.argv[2]
    output_germline = sys.argv[3]
    
    process_and_split_data(input_csv, output_somatic, output_germline)