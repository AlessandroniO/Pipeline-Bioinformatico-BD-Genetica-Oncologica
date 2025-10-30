import cbioportal_client
import pandas as pd
import argparse
import sys
import os # Importamos os para verificar la existencia de archivos

# ====================================================
# CONFIGURACIÓN DE ARGUMENTOS (PARA SNAKEMAKE)
# ====================================================
# Snakemake pasa la ruta del archivo de salida como primer argumento.
parser = argparse.ArgumentParser(description="Consulta variantes oncológicas frecuentes en cBioPortal y guarda el resultado en un CSV.")
parser.add_argument("output_file", help="Ruta del archivo CSV de salida especificado por Snakemake.")
args = parser.parse_args()

# ====================================================
# 1. Conexión y Configuración
# ====================================================
# ID del estudio Pan-Cancer (uno de los más grandes y completos)
study_id = "all_pan_cancer_atlas"
gene_panel_id = "msk_impact_500" 

try:
    # URL pública por defecto para la API
    client = cbioportal_client.ApiClient(api_url='https://www.cbioportal.org/api')
    print("Conexión exitosa a la API de cBioPortal.")
except Exception as e:
    print(f"Error al conectar a la API: {e}", file=sys.stderr)
    sys.exit(1) # Termina el script con un error si la conexión falla

# ====================================================
# 2. Hacer la Consulta
# ====================================================
try:
    # Obtener el resumen de mutaciones para el estudio y el panel de genes
    gene_summary = client.studies_get_mutated_genes(
        study_id=study_id,
        gene_panel_id=gene_panel_id,
        q_value_threshold=0.05 # Filtra por genes significativamente mutados
    )

    print(f"\nConsulta realizada al estudio: {study_id}")
    print(f"Total de genes significativamente mutados encontrados: {len(gene_summary)}")

    # ====================================================
    # 3. Procesar, Estructurar y Ordenar
    # ====================================================
    data = []
    for entry in gene_summary:
        data.append({
            "Gene": entry.gene.hugo_gene_symbol,
            "Frequency_Altered": entry.fraction_altered, 
            "Mutated_Cases": entry.mutated_cases_count,
            "Total_Cases": entry.all_cases_count,
            "Somatic_Mutation_Rate": entry.somatic_mutation_rate # Usado para ordenar
        })

    df = pd.DataFrame(data)
    # Ordenar por la Frecuencia de Mutación Somática
    df_sorted = df.sort_values(by="Somatic_Mutation_Rate", ascending=False)

    # ====================================================
    # 4. Guardar la Salida al Archivo (Paso CRÍTICO para Snakemake)
    # ====================================================
    # Guardamos el DataFrame completo en el archivo de salida especificado por Snakemake
    df_sorted.to_csv(args.output_file, index=False)
    print(f"\n✅ Resultados guardados exitosamente en: {args.output_file}")

    # (Opcional: puedes imprimir un resumen en consola si lo deseas, 
    # pero el archivo debe ser generado)
    print(df_sorted[['Gene', 'Somatic_Mutation_Rate', 'Mutated_Cases']].head(5).to_markdown(index=False))

except cbioportal_client.ApiException as e:
    print(f"Error en la consulta API: Asegúrate de que los IDs del estudio y del panel son correctos. {e}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"Ocurrió un error inesperado: {e}", file=sys.stderr)
    sys.exit(1)
