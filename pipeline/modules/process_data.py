import sys
import pandas as pd
import os

def process_data():
    if len(sys.argv) < 3:
        print("ERROR: Faltan argumentos. Se requiere archivo de entrada y salida.")
        sys.exit(1)
        
    input_csv = sys.argv[1] 
    output_csv = sys.argv[2] 

    print(f"INFO: Leyendo datos de: {input_csv}")
    
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        print(f"❌ ERROR al leer el archivo CSV: {e}")
        sys.exit(1)

    initial_count = len(df)
    
    # =========================================================================
    # LÓGICA DE LIMPIEZA Y CORRECCIÓN DE MAPEO (CORREGIDO)
    # =========================================================================

    # 1. Renombrar Frec_alelica y la columna que contiene la Significancia Clínica (Detalle_variante)
    df = df.rename(columns={'Detalle_variante': 'Significancia_clinica_CORREGIDA',
                            'Frec_alelica': 'AF_Alelica'})

    # 2. Limpieza Crítica: Eliminar variantes que no tienen Muestra ni Tumor asociado.
    df_filtered = df.dropna(subset=['Id_muestra', 'Id_tipo_tumor'], how='all')

    # 3. Imputación: Si la significancia (corregida) es nula (por ejemplo, en las primeras filas de variantes_extraidas.csv), imputamos "Significancia_Desconocida".
    df_final = df_filtered.copy()
    df_final['Significancia_clinica_CORREGIDA'] = df_final['Significancia_clinica_CORREGIDA'].fillna('Significancia_Desconocida')
    
    # 4. Asegurar que el nombre de columna final coincida con el script de VCF
    df_final = df_final.rename(columns={'Significancia_clinica_CORREGIDA': 'Significancia_clinica'})
    
    # 5. Eliminar columnas ID que no se usarán para el análisis final
    cols_to_drop = [col for col in df_final.columns if col.startswith('Id_') and col != 'Id_muestra']
    df_final = df_final.drop(columns=cols_to_drop, errors='ignore')
    
    final_count = len(df_final)
    
    print(f"INFO: Filas iniciales: {initial_count}")
    print(f"INFO: Filas descartadas por falta de info crítica: {initial_count - final_count}")
    print(f"INFO: Total de filas limpias y válidas: {final_count}")
    
    # Exportar el archivo limpio
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df_final.to_csv(output_csv, index=False)
    print(f"✅ Éxito: Datos limpios exportados a {output_csv}")


if __name__ == "__main__":
    process_data()