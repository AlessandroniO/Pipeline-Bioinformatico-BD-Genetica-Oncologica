import mysql.connector
import pandas as pd
import sys
import os
from cyvcf2 import VCF
from dotenv import load_dotenv

# =========================================================================
# 1. CONFIGURACIÓN Y CONEXIÓN
# =========================================================================

load_dotenv(os.path.join(os.path.dirname(__file__), '..', '..', '.env'))

DB_CONFIG = {
    'host': os.getenv("MYSQL_HOST"),
    'user': os.getenv("MYSQL_USER"),
    'password': os.getenv("MYSQL_ROOT_PASSWORD"),
    'database': "dw_banco_de_datos", 
    'port': os.getenv("MYSQL_PORT", 3306)
}

BATCH_SIZE = 1000 
TABLE_NAME = "Variantes"        # Tabla final
STAGING_TABLE_NAME = "Clinvar_Staging" # 💥 Nueva tabla intermedia

def get_vcf_path_and_output_flag():
    if len(sys.argv) < 3:
        print("❌ WARNING: Usando rutas de fallback para pruebas.")
        return "GRCh38_latest_clinvar.vcf", "data/intermediate/db_enriched.flag"
    return sys.argv[1], sys.argv[2]

# =========================================================================
# 2. FUNCIÓN DE INSERCIÓN POR LOTES (Carga a la tabla de Staging)
# =========================================================================

def insert_data_in_batches(conn, cursor, df, table_name):
    """Inserta el DataFrame en la tabla de staging usando ON DUPLICATE KEY UPDATE."""
    
    # Mapeo del DataFrame a las columnas de la tabla de staging
    df_columns_to_map = [
        "chrom", "pos", "ref", "alt", "qual", "filter", "tipo_variante", "source", "id"
    ]
    
    # Columnas correspondientes en la tabla SQL 'Clinvar_Staging'
    sql_columns = [
        "Cromosoma", "Posicion", "Alelo_referencia", "Alelo_alternativo", 
        "Calidad", "Filtro", "Tipo_variante", "Fuente", "Variante"
    ]
    
    placeholders = ', '.join(['%s'] * len(sql_columns))
    
    # 💥 Usa ON DUPLICATE KEY UPDATE (Requiere el índice UNIQUE en Clinvar_Staging)
    insert_query = f"""
        INSERT INTO {table_name} ({', '.join(sql_columns)}) 
        VALUES ({placeholders})
        ON DUPLICATE KEY UPDATE
        Calidad = VALUES(Calidad),
        Filtro = VALUES(Filtro),
        Tipo_variante = VALUES(Tipo_variante),
        Fuente = VALUES(Fuente),
        Variante = VALUES(Variante)
    """

    print(f"INFO: Cargando {len(df)} registros en '{table_name}' en lotes de {BATCH_SIZE}...")
    total_processed = 0
    df_insert = df[df_columns_to_map]

    for i in range(0, len(df_insert), BATCH_SIZE):
        batch_df = df_insert.iloc[i:i + BATCH_SIZE]
        batch_data = [tuple(x) for x in batch_df.to_numpy()]
        
        try:
            cursor.executemany(insert_query, batch_data)
            conn.commit()
            total_processed += len(batch_df)
        except mysql.connector.Error as err:
            print(f"❌ ERROR al insertar lote {i // BATCH_SIZE + 1}: {err}")
            conn.rollback() 
            raise 
            
    print(f"✅ Carga en {table_name} completada. Total procesado: {total_processed}")

# =========================================================
# 3. FUNCIÓN DE FUSIÓN (Actualiza e Inserta en Variantes)
# =========================================================

def merge_staging_to_variantes(conn, cursor, staging_table, final_table, batch_size=50000):
    """Fusión incremental entre Clinvar_Staging y Variantes."""
    print(f"\nINFO: Iniciando fusión incremental de {staging_table} a {final_table}...")

    # Obtener todas las variantes de staging en memoria ligera
    cursor.execute(f"SELECT Variante FROM {staging_table}")
    all_variants = [row[0] for row in cursor.fetchall()]
    total_variants = len(all_variants)
    print(f"INFO: Total de variantes en staging: {total_variants}")

    updated_total = 0
    inserted_total = 0

    for start in range(0, total_variants, batch_size):
        batch = all_variants[start:start + batch_size]
        placeholders = ', '.join(['%s'] * len(batch))
        
        print(f"🧩 Procesando batch {start // batch_size + 1} ({len(batch)} variantes)...")

        # === 1. UPDATE parcial ===
        update_query = f"""
            UPDATE {final_table} AS v
            INNER JOIN {staging_table} AS cs 
                ON v.Variante = cs.Variante
            SET
                v.Cromosoma = cs.Cromosoma,
                v.Posicion = cs.Posicion,
                v.Alelo_referencia = cs.Alelo_referencia,
                v.Alelo_alternativo = cs.Alelo_alternativo,
                v.Calidad = cs.Calidad,
                v.Filtro = cs.Filtro,
                v.Tipo_variante = cs.Tipo_variante,
                v.Fuente = cs.Fuente
            WHERE v.Cromosoma IS NULL
              AND cs.Variante IN ({placeholders});
        """
        cursor.execute(update_query, batch)
        updated_total += cursor.rowcount
        conn.commit()

        # === 2. INSERT parcial ===
        insert_query = f"""
            INSERT INTO {final_table} (
                Variante, Tipo_modificacion, Descripcion,
                Cromosoma, Posicion, Alelo_referencia, Alelo_alternativo, 
                Calidad, Filtro, Tipo_variante, Fuente
            )
            SELECT 
                cs.Variante, 
                'single nucleotide variant', 
                'Data loaded from ClinVar VCF',
                cs.Cromosoma, cs.Posicion, cs.Alelo_referencia, cs.Alelo_alternativo,
                cs.Calidad, cs.Filtro, cs.Tipo_variante, cs.Fuente
            FROM {staging_table} AS cs
            LEFT JOIN {final_table} AS v ON cs.Variante = v.Variante
            WHERE v.Id_variante IS NULL
              AND cs.Variante IN ({placeholders});
        """
        cursor.execute(insert_query, batch)
        inserted_total += cursor.rowcount
        conn.commit()

        print(f"✅ Batch {start // batch_size + 1} -> {cursor.rowcount} filas insertadas / {updated_total} actualizadas")

    print(f"🏁 Fusión completada: {updated_total} filas actualizadas, {inserted_total} filas insertadas en total.")


# =========================================================================
# 4. FUNCIÓN PRINCIPAL run_enrichment (FLUJO ETL COMPLETO)
# =========================================================================

def run_enrichment():
    VCF_PATH, OUTPUT_FLAG = get_vcf_path_and_output_flag()

    print(f"INFO: Conectando a la BD '{DB_CONFIG['database']}'...")
    try:
        conn = mysql.connector.connect(**DB_CONFIG)
        cursor = conn.cursor()
    except mysql.connector.Error as err:
        print(f"❌ ERROR: Fallo de conexión a MySQL: {err}")
        sys.exit(1)

    # ... (CÓDIGO DE LECTURA DE VCF) ...
    if not os.path.exists(VCF_PATH):
        print(f"❌ ERROR: No se encuentra el archivo ClinVar VCF en {VCF_PATH}")
        sys.exit(1)

    print(f"INFO: Leyendo datos desde ClinVar: {VCF_PATH}")
    try:
        vcf = VCF(VCF_PATH) 
    except Exception as e:
        print(f"❌ ERROR al leer VCF: {e}")
        conn.close()
        sys.exit(1)

    records = []
    for i, variant in enumerate(vcf):
        chrom_val = str(variant.CHROM).replace('chr', '') if variant.CHROM else "NA"
        try:
            pos_val = int(variant.POS)
        except (TypeError, ValueError):
            pos_val = 0
        qual_val = float(variant.QUAL) if variant.QUAL is not None else 0.0

        records.append({
            "chrom": chrom_val,
            "pos": pos_val,
            "id": variant.ID or f"CV_{i}", 
            "ref": variant.REF,
            "alt": ",".join(variant.ALT),
            "qual": qual_val,
            "filter": ";".join(variant.FILTERS) if variant.FILTERS else "PASS",
            "tipo_variante": variant.INFO.get("CLNVC", "ND"),
            "source": "ClinVar" 
        })
    
    df = pd.DataFrame(records)

    # Limpieza consistente 
    df["ref"] = df["ref"].fillna("N")
    df["alt"] = df["alt"].fillna("N")
    df["qual"] = df["qual"].fillna(0.0)
    df["filter"] = df["filter"].fillna("UNKNOWN")
    df["tipo_variante"] = df["tipo_variante"].fillna("ND")
    
    # 💥 FLUJO ETL 💥
    
    # 1. Truncar la tabla de staging
    try:
        print(f"INFO: Truncando tabla {STAGING_TABLE_NAME} para limpieza...")
        cursor.execute(f"TRUNCATE TABLE {STAGING_TABLE_NAME}")
        conn.commit()
    except mysql.connector.Error as err:
        # Esto puede fallar si la tabla aún no existe, lo cual es manejable si el Snakefile la crea primero.
        print(f"❌ WARNING: Fallo al truncar {STAGING_TABLE_NAME} (podría no existir): {err}")
        # NO salir, ya que la inserción fallaría si la tabla no existe.

    # 2. Cargar los datos del VCF a la tabla STAGING
    try:
        insert_data_in_batches(conn, cursor, df, STAGING_TABLE_NAME)
    except Exception as e:
        print(f"❌ ERROR FATAL durante la carga de staging: {e}")
        conn.close()
        sys.exit(1)

    # 3. Fusionar datos de STAGING a la tabla FINAL Variantes
    try:
        merge_staging_to_variantes(conn, cursor, STAGING_TABLE_NAME, TABLE_NAME)
    except Exception as e:
        print(f"❌ ERROR FATAL durante la fusión de datos: {e}")
        conn.rollback()
        conn.close()
        sys.exit(1)

    # =========================================================
    # EXPORTAR CSV INTERMEDIO (código sin cambios)
    # =========================================================
    try:
        output_csv_path = "data/intermediate/vcf_pure_clinvar.csv"
        os.makedirs(os.dirname(output_csv_path), exist_ok=True)
        df.to_csv(output_csv_path, index=False)
        print(f"✅ Archivo intermedio generado correctamente: {output_csv_path}")
    except Exception as e:
        print(f"❌ ERROR al generar CSV intermedio: {e}")
        sys.exit(1)

    # =========================================================
    # CREAR BANDERA DE ÉXITO (código sin cambios)
    # =========================================================
    try:
        os.makedirs(os.path.dirname(OUTPUT_FLAG), exist_ok=True)
        with open(OUTPUT_FLAG, "w") as f:
            f.write("Enriquecimiento ClinVar completado.\n")
        print(f"✅ Bandera de éxito creada en: {OUTPUT_FLAG}")
    except Exception as e:
        print(f"❌ ERROR al crear bandera de salida: {e}")
        sys.exit(1)

    cursor.close()
    conn.close()
    print("✅ Proceso finalizado correctamente.")

# =========================================================================
# 5. EJECUCIÓN (código sin cambios)
# =========================================================================

if __name__ == "__main__":
    run_enrichment()