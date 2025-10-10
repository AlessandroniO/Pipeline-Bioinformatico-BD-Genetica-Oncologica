import sys
import pandas as pd
import os
import re

# =========================================================================
# 1. FUNCIONES AUXILIARES
# =========================================================================

def infer_genomic_type(tipo_muestra):
    """
    Infiere el tipo de variante genómica (SOMATIC/GERMLINE) basándose
    en la fuente física de la muestra proporcionada en el CSV de origen.
    """
    if pd.isna(tipo_muestra):
        return 'UNKNOWN'
    
    # Limpieza y estandarización del texto para la inferencia
    tipo_muestra_str = str(tipo_muestra).strip().upper()

    # Si la muestra es Sangre (Blood), se asume que es la muestra de control.
    if 'SANGRE' in tipo_muestra_str or 'BLOOD' in tipo_muestra_str:
        return 'GERMLINE'
    # Cualquier otra fuente (tejido, líquido, etc.) se asume tumoral.
    elif 'TEJIDO' in tipo_muestra_str or 'LIQUIDO' in tipo_muestra_str or 'FLUID' in tipo_muestra_str or 'TISSUE' in tipo_muestra_str:
        return 'SOMATIC'
    else:
        return 'UNKNOWN'

# =========================================================================
# 2. FUNCIÓN PRINCIPAL DE PROCESAMIENTO
# =========================================================================

def process_and_split_data(input_csv, output_somatic_csv, output_germline_csv):
    """
    Carga el CSV, realiza la limpieza, estandarización e inferencia del tipo
    genómico, y divide los datos en archivos somáticos y germinales.
    """
    try:
        # 1. CARGA DE DATOS
        df = pd.read_csv(input_csv)
        print(f"INFO: Datos cargados. Filas iniciales: {len(df)}")

        # 2. LIMPIEZA Y ESTANDARIZACIÓN DE DATOS
        
        # Rellenar valores nulos/vacíos en metadatos clave con 'DESCONOCIDO' o 0
        # CRUCIAL: Usamos los nombres de columna cortos ('tipo_tumor', 'sexo', 'barrio', 'edad')
        # ya que db_extraction.py se encarga de renombrarlos con alias.
        
        # Limpieza de metadatos del paciente
        df['tipo_tumor'] = df['tipo_tumor'].fillna('DESCONOCIDO')
        df['sexo'] = df['sexo'].fillna('DESCONOCIDO')
        
        # Corregir la edad: Si es nulo, poner 0 y convertir a entero.
        df['edad'] = df['edad'].fillna(0).astype(int)
        
        # Corregir el barrio:
        # ERROR CORREGIDO: Ahora usamos 'barrio' directamente, no el nombre largo.
        df['barrio'] = df['barrio'].fillna('DESCONOCIDO') 

        # Estandarización de la columna 'chrom'
        df['chrom'] = df['chrom'].astype(str).str.replace('^chr', '', regex=True)
        df['chrom'] = 'chr' + df['chrom']

        # 3. INFERENCIA DEL TIPO GENÓMICO (SOMATIC/GERMLINE)
        
        # Creamos una nueva columna 'tipo_genomico' basada en la columna 'tipo_de_muestra'
        df['tipo_genomico'] = df['tipo_de_muestra'].apply(infer_genomic_type)
        
        # Verificación de la división (opcional, pero útil)
        somatic_count = (df['tipo_genomico'] == 'SOMATIC').sum()
        germline_count = (df['tipo_genomico'] == 'GERMLINE').sum()
        unknown_count = (df['tipo_genomico'] == 'UNKNOWN').sum()

        print(f"INFO: Inferencia completada. Somáticas: {somatic_count}, Germinales: {germline_count}, Desconocidas: {unknown_count}")

        # 4. DIVISIÓN Y EXPORTACIÓN
        
        # Creamos los sub-DataFrames para somático y germinal
        df_somatic = df[df['tipo_genomico'] == 'SOMATIC']
        df_germline = df[df['tipo_genomico'] == 'GERMLINE']
        
        # Crear directorios de salida si no existen
        os.makedirs(os.path.dirname(output_somatic_csv), exist_ok=True)
        
        # Exportar los archivos limpios
        df_somatic.to_csv(output_somatic_csv, index=False)
        print(f"INFO: ✅ Datos Somáticos exportados a {output_somatic_csv} ({len(df_somatic)} filas)")
        
        df_germline.to_csv(output_germline_csv, index=False)
        print(f"INFO: ✅ Datos Germinales exportados a {output_germline_csv} ({len(df_germline)} filas)")

    except FileNotFoundError:
        print(f"❌ ERROR: El archivo de entrada no se encontró: {input_csv}")
        sys.exit(1)
    except KeyError as e:
        print(f"❌ ERROR CRÍTICO: Columna no encontrada. El archivo de entrada debe tener el formato correcto.")
        print(f"COLUMNA FALTANTE: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error inesperado durante el procesamiento de datos: {e}")
        sys.exit(1)


# =========================================================================
# 3. EJECUCIÓN DEL SCRIPT
# =========================================================================

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python data_processing.py <input_csv> <output_somatic_csv> <output_germline_csv>")
        sys.exit(1)
        
    input_csv = sys.argv[1]
    output_somatic_csv = sys.argv[2]
    output_germline_csv = sys.argv[3]
    
    process_and_split_data(input_csv, output_somatic_csv, output_germline_csv)
