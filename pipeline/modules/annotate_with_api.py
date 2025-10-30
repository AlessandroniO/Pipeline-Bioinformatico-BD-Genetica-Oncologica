import pandas as pd
import requests
import sys
import time
import re
import os # Necesario para la comprobación final de escritura

# Endpoint de Mutalyzer para normalizar HGVS
MUTALYZER_API_URL = "https://mutalyzer.nl/api/normalize"

def get_vcf_coordinates(hgvs_variant):
    """
    Consulta la API de Mutalyzer para obtener las coordenadas VCF (HG38)
    a partir de la notación HGVS.
    """
    
    # 1. Limpieza y preparación de datos: Mutalyzer necesita solo el ID y la variante c.
    # Ejemplo: NM_001162498.3:c.105T>C
    clean_hgvs = hgvs_variant.split(' ')[0]

    # 2. Configuración de la consulta
    headers = {"Content-Type": "application/json"}
    
    try:
        time.sleep(0.1) 
        
        # Mutalyzer toma una lista de descripciones
        response = requests.post(MUTALYZER_API_URL, 
                                 json={"descriptions": [clean_hgvs]}, 
                                 headers=headers, 
                                 timeout=20)
        response.raise_for_status()

        data = response.json()
        
        # 3. Procesamiento de la respuesta de Mutalyzer
        if not data or 'error' in data[0] or 'gDNA' not in data[0]:
            # Mutalyzer falló o no encontró mapeo gDNA
            return {'#CHROM': 'NA_MUT', 'POS': 'NA_MUT', 'REF': 'NA_MUT', 'ALT': 'NA_MUT'}
            
        result = data[0]
        g_notation = result.get('gDNA') # Ej: 'NC_000017.11:g.7668407G>C'
        
        # 4. Extracción de coordenadas y alelos usando expresiones regulares (BLOQUE CRÍTICO)
        try:
            # Búsqueda más simple: cualquier número de cromosoma y POS, REF, ALT
            # Esto captura cualquier notación g. que contenga POS, REF, y ALT
            chrom_match = re.search(r'(NC_0000(\d+))', g_notation) # Captura el cromosoma
            vcf_match = re.search(r'g\.(\d+)([A-Z])>([A-Z])', g_notation) # Captura POS, REF, ALT
            
            if chrom_match and vcf_match:
                # Obtenemos solo el número de cromosoma y lo limpiamos (Ej: '13')
                chrom_num = chrom_match.group(2).lstrip('0')
                
                # Extraer POS, REF, ALT
                pos, ref, alt = vcf_match.groups()
                
                return {
                    '#CHROM': f'chr{chrom_num}', # Añadimos 'chr' para el formato VCF
                    'POS': int(pos),
                    'REF': ref,
                    'ALT': alt
                }
            
            # Si el regex no encuentra coincidencias, se clasifica como NA
            return {'#CHROM': 'NA_MUT', 'POS': 'NA_MUT', 'REF': 'NA_MUT', 'ALT': 'NA_MUT'}
            
        except Exception as e:
            # Captura cualquier error en el análisis de la cadena gDNA (fallo de código)
            print(f"Error de parsing interno para gDNA: {g_notation}. Error: {e}", file=sys.stderr)
            return {'#CHROM': 'PARSE_FAIL', 'POS': 'PARSE_FAIL', 'REF': 'PARSE_FAIL', 'ALT': 'PARSE_FAIL'}

    except requests.exceptions.RequestException as e:
        print(f"Error de red o API para variante {clean_hgvs}: {e}", file=sys.stderr)
    except Exception as e:
        print(f"Error inesperado al procesar la respuesta para {clean_hgvs}: {e}", file=sys.stderr)
        
    return {'#CHROM': 'FAIL', 'POS': 'FAIL', 'REF': 'FAIL', 'ALT': 'FAIL'}


def annotate_csv(input_csv, output_csv):
    """Procesa el CSV, consulta la API y guarda el resultado."""
    print(f"Iniciando anotación por API (Mutalyzer) para {input_csv}...")
    
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        print(f"Error al leer el archivo CSV: {e}", file=sys.stderr)
        sys.exit(1)

    # 1. Identificar variantes únicas para optimizar llamadas API
    unique_variants = df['variante'].unique()
    
    # 2. Realizar las consultas API para obtener el mapeo VCF
    variant_map = {}
    num_variants = len(unique_variants)
    
    for i, hgvs in enumerate(unique_variants):
        if hgvs not in variant_map:
            vcf_data = get_vcf_coordinates(hgvs)
            variant_map[hgvs] = vcf_data
            
            if (i + 1) % 10 == 0 or (i + 1) == num_variants:
                 print(f"-> Procesadas {i + 1}/{num_variants} variantes.")
                 
    
    print(f"Anotación de {num_variants} variantes completada. Fusionando datos.")

    # 3. Preparar el DataFrame con las coordenadas
    coord_cols = ['#CHROM', 'POS', 'REF', 'ALT']
    
    # Crear una columna de mapeo con la HGVS limpia (solo la notación de transcripto)
    df['hgvs_clean'] = df['variante'].str.split(' ').str[0]
    
    # Crear un DataFrame de mapeo a partir del diccionario de resultados
    map_df = pd.DataFrame(variant_map).T.reset_index(names=['hgvs_clean'])
    
    # 4. Fusionar los datos VCF con el DataFrame original
    df_merged = df.merge(map_df, on='hgvs_clean', how='left')
    
    # Limpieza final
    df_merged = df_merged.drop(columns=['hgvs_clean'])

    print(f"Éxito. {df_merged.shape[0]} filas anotadas. Guardando en {output_csv}")
    
    # 5. Guardar el archivo con manejo de errores de escritura
    try:
        df_merged.to_csv(output_csv, index=False)
    except Exception as e:
        print(f"ERROR CRÍTICO: No se pudo guardar el archivo {output_csv}. Revise permisos de Docker. Error: {e}", file=sys.stderr)
        # Forzamos la salida con error para que Snakemake falle correctamente
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python annotate_with_api.py <input_csv> <output_csv>", file=sys.stderr)
        sys.exit(1)

    input_csv = sys.argv[1]
    output_csv = sys.argv[2]
    
    # Se añade un try/except general para capturar cualquier error de nivel superior
    try:
        annotate_csv(input_csv, output_csv)
    except Exception as e:
        print(f"ERROR CRÍTICO: Fallo en la función principal annotate_csv: {e}", file=sys.stderr)
        sys.exit(1)