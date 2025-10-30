import sys
import os
import pandas as pd
# Aseguramos que cyvcf2 esté disponible, ya que es mucho más rápido que pysam o pyvcf
try:
    from cyvcf2 import VCF
except ImportError:
    print("❌ ERROR: La librería 'cyvcf2' no está instalada. Ejecute 'pip install cyvcf2' o revise su entorno.")
    sys.exit(1)

# =========================================================================
# 1. FUNCIÓN AUXILIAR DE RUTA
# =========================================================================
def get_io_paths():
    """
    Obtiene las rutas de entrada (VCF) y salida (CSV) de los argumentos
    de línea de comandos (proporcionados por Snakemake).
    """
    if len(sys.argv) != 3:
        print("Uso: python extract_vcf_data.py <input_vcf_path> <output_csv_path>")
        # En el contexto de Snakemake, la ruta VCF será 'GRCh38_latest_clinvar.vcf.gz'
        return "GRCh38_latest_clinvar.vcf.gz", "data/intermediate/vcf_only_extracted.csv"
    
    return sys.argv[1], sys.argv[2]

# =========================================================================
# 2. FUNCIÓN PRINCIPAL DE EXTRACCIÓN VCF
# =========================================================================
def extract_vcf_data():
    """
    Lee un archivo VCF, extrae las columnas fundamentales (CHROM, POS, ID, REF, ALT, QUAL, FILTER)
    y lo guarda en un CSV intermedio.
    """
    INPUT_VCF, OUTPUT_CSV = get_io_paths()
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    
    print(f"INFO: Leyendo datos de ClinVar desde: {INPUT_VCF}")

    if not os.path.exists(INPUT_VCF):
        print(f"❌ ERROR CRÍTICO: Archivo VCF no encontrado en {INPUT_VCF}")
        sys.exit(1)

    records = []
    try:
        vcf = VCF(INPUT_VCF)

        for i, variant in enumerate(vcf):
            # Usamos una clave compuesta (CHROM, POS, REF, ALT) como un identificador robusto.
            # Además, incluimos el ID de ClinVar si existe.
            records.append({
                "chrom": variant.CHROM,
                "pos": variant.POS,
                "ref": variant.REF,
                "alt": ",".join(variant.ALT),
                "id_clinvar": variant.ID or f"CV_{i}", # Usar ID o un ID generado si no tiene
                "qual": variant.QUAL if variant.QUAL is not None else 0.0,
                "filter": ";".join(variant.FILTERS) if variant.FILTERS else "PASS",
                # Opcional: Extraer la clase de variante clínica (CLNVC) si existe
                "clnvc_type": variant.INFO.get("CLNVC", None)
            })
        
        df = pd.DataFrame(records)
        print(f"INFO: Se extrajeron {len(df)} variantes VCF puras.")

        # Guardar el DataFrame
        df.to_csv(OUTPUT_CSV, index=False)
        print(f"INFO: ✅ Datos VCF exportados a {OUTPUT_CSV}")

    except Exception as e:
        print(f"❌ ERROR inesperado durante la lectura o procesamiento del VCF: {e}")
        sys.exit(1)

# =========================================================================
# 3. EJECUCIÓN
# =========================================================================
if __name__ == "__main__":
    extract_vcf_data()
