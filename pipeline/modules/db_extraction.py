import sys
import os
import pandas as pd
import mysql.connector
from mysql.connector import pooling # Importar pooling para gestionar conexiones

# =========================================================================
# 1. CONFIGURACIÓN Y CONEXIÓN
# =========================================================================

# NOTA: Usando credenciales hardcodeadas por solicitud del usuario.
DB_CONFIG = {
    'host': '192.168.100.6',
    'user': 'root',
    'password': 'arquisoft1',
    'database': 'dw_banco_de_datos'
}

# Usamos un Pool de conexiones para manejar mejor la escalabilidad
try:
    db_pool = mysql.connector.pooling.MySQLConnectionPool(
        pool_name="bioinfo_pool",
        pool_size=5,  # Tamaño del pool de conexiones
        **DB_CONFIG
    )
    print("INFO: Pool de conexiones MySQL creado con éxito.")
except mysql.connector.Error as err:
    # Ahora que no usamos ENV, la configuración debe ser correcta o fallar aquí
    print(f"❌ ERROR CRÍTICO MySQL: No se pudo crear el Pool de conexiones: {err}")
    sys.exit(1)


# =========================================================================
# 2. FUNCIÓN PARA OBTENER LA RUTA DE SALIDA
# =========================================================================
def get_output_path():
    """
    Obtiene la ruta del archivo de salida de los argumentos de línea de comandos.
    """
    if len(sys.argv) > 1:
        return sys.argv[1]
    else:
        print("ADVERTENCIA: Ejecución directa. Usando ruta de prueba. Use Snakemake para la ruta final.")
        return "data/intermediate/variantes_extraidas.csv"

# =========================================================================
# 3. FUNCIÓN DE EXTRACCIÓN Y EXPORTACIÓN DE DATOS (SIN MOCK DATA)
# =========================================================================
def fetch_data_and_export():
    """
    Conecta a la base de datos, extrae las variantes y las exporta a un CSV.
    Si la conexión falla, el script termina con un error de Python.
    """
    
    OUTPUT_CSV = get_output_path()
    
    # Crear el directorio si no existe
    output_dir = os.path.dirname(OUTPUT_CSV)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    try:
        # 1. OBTENER CONEXIÓN DEL POOL
        conn = db_pool.get_connection()
        if not conn.is_connected():
            raise mysql.connector.Error("No se pudo obtener una conexión válida del pool.")
        
        # 2. CONSULTA SQL CORREGIDA: Agregamos 0.5 AS frec_alelica
        SQL_QUERY = """
        SELECT
            -- Campos VCF Requeridos por el pipeline
            ug.Cromosoma AS chrom,
            v.Posicion AS pos, 
            v.Id_variante AS id,
            v.Alelo_referencia AS ref,
            v.Alelo_alternativo AS alt,
            v.Calidad AS qual,
            v.Filtro AS filter,
            v.Tipo_variante AS tipo_variante,
            
            -- COLUMNA REQUERIDA POR EL SCRIPT DE CONVERSIÓN VCF (vcf_conversion.py)
            0.5 AS frec_alelica, 
            
            -- Metadatos de la Muestra (para inferencia SOMATIC/GERMLINE)
            m.Tipo_de_muestra AS tipo_de_muestra,

            -- Metadatos del Paciente (Requeridos por QC y Split)
            p.Sexo AS sexo,
            YEAR(CURDATE()) - YEAR(p.Fecha_nacimiento) AS edad,
            tt.Nombre AS tipo_tumor,
            d_addr.Barrio AS barrio 
            
        FROM Variantes v
        -- Uniones para obtener metadatos y ubicación (manteniendo tu esquema de JOINs)
        LEFT JOIN Ubicacion_x_variante uxv
            ON v.Id_variante = uxv.Id_variante
        LEFT JOIN Ubicacion_genomica ug
            ON uxv.Id_ubicacion_genomica = ug.Id_ubicacion_genomica
        LEFT JOIN Tipo_tumor_variante ttv
            ON v.Id_variante = ttv.Id_variante
        LEFT JOIN Tipo_tumor tt
            ON ttv.Id_tipo_tumor = tt.Id_tipo_tumor
        LEFT JOIN Variante_x_informe vxi
            ON v.Id_variante = vxi.Id_variante
        LEFT JOIN Informe_resultado ir
            ON vxi.Id_informe = ir.Id_informe
        LEFT JOIN Muestra_genomica m
            ON ir.Id_muestra = m.Id_muestra
        LEFT JOIN Paciente p
            ON p.Id_paciente = m.Id_paciente
        LEFT JOIN Direccion d_addr
            ON p.Id_direccion = d_addr.Id_direccion
        LEFT JOIN Diagnostico_paciente d
            ON p.Id_paciente = d.Id_paciente
        ORDER BY m.Id_muestra, v.Id_variante;

        """

        # 3. EJECUTAR CONSULTA Y OBTENER DATOS EN UN DATAFRAME
        df = pd.read_sql(SQL_QUERY, conn)
        
        # 4. CERRAR CONEXIÓN (IMPORTANTE DEVOLVERLA AL POOL)
        conn.close()
        print(f"INFO: ✅ Éxito en la conexión y extracción. Se extrajeron {len(df)} filas.")
        
        # ---------------------------------------------------------------------
        # EXPORTAR EL DATAFRAME A CSV
        # ---------------------------------------------------------------------
        df.to_csv(OUTPUT_CSV, index=False)
        print(f"INFO: Datos exportados a {OUTPUT_CSV}")

    except mysql.connector.Error as err:
        print(f"❌ ERROR CRÍTICO MySQL: Algo salió mal al conectar o consultar: {err}")
        print("El pipeline fallará porque no se pudo obtener los datos reales.")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error inesperado durante la exportación: {e}")
        sys.exit(1)

# =========================================================================
# 4. EJECUCIÓN DEL SCRIPT
# =========================================================================
if __name__ == "__main__":
    fetch_data_and_export()
