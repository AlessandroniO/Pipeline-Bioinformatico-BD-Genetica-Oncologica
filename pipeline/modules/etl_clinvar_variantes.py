import mysql.connector
import pandas as pd
import sys
import os
from cyvcf2 import VCF

DB_CONFIG = {
    'host': '192.168.100.6',
    'user': 'root',
    'password': 'arquisoft1',
    'database': 'dw_banco_de_datos'
}

VCF_PATH = "/data/GRCh38_latest_clinvar.vcf.gz"

def run_enrichment():
    print(f"INFO: Conectando a la BD '{DB_CONFIG['database']}' para enriquecimiento...")
    try:
        conn = mysql.connector.connect(**DB_CONFIG)
        cursor = conn.cursor()
    except mysql.connector.Error as err:
        print(f"❌ ERROR: Fallo de conexión a MySQL: {err}")
        sys.exit(1)

    if not os.path.exists(VCF_PATH):
        print(f"❌ ERROR: No se encuentra el archivo ClinVar VCF en {VCF_PATH}")
        sys.exit(1)

    print(f"INFO: Leyendo datos desde ClinVar: {VCF_PATH}")
    vcf = VCF(VCF_PATH)
    records = []

    for i, variant in enumerate(vcf):
        if i >= 5000:  # evitar saturar MySQL
            break
        records.append({
            "chrom": variant.CHROM,
            "pos": variant.POS,
            "id": variant.ID or f"CV_{i}",
            "ref": variant.REF,
            "alt": ",".join(variant.ALT),
            "qual": float(variant.QUAL) if variant.QUAL is not None else 0.0,
            "filter": ";".join(variant.FILTERS) if variant.FILTERS else "PASS",
            "type": variant.INFO.get("CLNVC", "SNV")
        })
    df = pd.DataFrame(records)
    print(f"INFO: Se leyeron {len(df)} variantes de ClinVar.")

    # --- ACTUALIZAR variantes existentes ---
    print(f"INFO: Actualizando variantes existentes...")
    update_query = """
        UPDATE Variantes
        SET Cromosoma=%s, Posicion=%s, Alelo_referencia=%s,
            Alelo_alternativo=%s, Calidad=%s, Filtro=%s, Tipo_variante=%s
        WHERE Variante = %s OR Tipo_modificacion = %s
    """
    update_count = 0
    for _, row in df.iterrows():
        cursor.execute(update_query, (
            row["chrom"], row["pos"], row["ref"], row["alt"],
            row["qual"], row["filter"], row["type"], row["id"], row["id"]
        ))
        update_count += cursor.rowcount
    conn.commit()
    print(f"INFO: ✅ {update_count} variantes actualizadas.")

    # --- INSERTAR variantes nuevas ---
    print("INFO: Insertando variantes nuevas desde ClinVar...")
    insert_query = """
        INSERT INTO Variantes (Variante, Tipo_modificacion, Cromosoma, Posicion,
                               Alelo_referencia, Alelo_alternativo, Calidad, Filtro, Tipo_variante, Fuente)
        SELECT %s, %s, %s, %s, %s, %s, %s, %s, %s, 'ClinVar'
        WHERE NOT EXISTS (
            SELECT 1 FROM Variantes 
            WHERE Cromosoma=%s AND Posicion=%s AND Alelo_referencia=%s AND Alelo_alternativo=%s
        )
    """
    insert_count = 0
    for _, row in df.iterrows():
        cursor.execute(insert_query, (
            row["id"], row["id"], row["chrom"], row["pos"], row["ref"],
            row["alt"], row["qual"], row["filter"], row["type"],
            row["chrom"], row["pos"], row["ref"], row["alt"]
        ))
        insert_count += cursor.rowcount
    conn.commit()
    print(f"INFO: ✅ {insert_count} nuevas variantes insertadas desde ClinVar.")
    print("✅ Proceso de Enriquecimiento completado.")

    cursor.close()
    conn.close()

if __name__ == "__main__":
    run_enrichment()
