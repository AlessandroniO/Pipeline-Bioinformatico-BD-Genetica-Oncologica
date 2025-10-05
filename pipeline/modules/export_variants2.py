import sys
import os
import pandas as pd
import mysql.connector

# =========================================================================
# 1. FUNCIÓN PARA OBTENER LA RUTA DE SALIDA
# =========================================================================
def get_output_path():
    """
    Obtiene la ruta del archivo de salida de los argumentos de línea de comandos.
    """
    if len(sys.argv) > 1:
        return sys.argv[1]
    else:
        # Modo de prueba si se ejecuta sin Snakemake
        print("ADVERTENCIA: Ejecución directa. Usando ruta de prueba. Use Snakemake para la ruta final.")
        return "data/processed/variantes_temp.csv"

# =========================================================================
# 2. FUNCIÓN DE EXTRACCIÓN Y EXPORTACIÓN DE DATOS
# =========================================================================
def fetch_data_and_export():
    """
    Conecta a la base de datos, extrae las variantes y las exporta a un CSV.
    """
    
    # Obtener la ruta de salida.
    OUTPUT_CSV = get_output_path()
    
    # Crear el directorio si no existe
    output_dir = os.path.dirname(OUTPUT_CSV)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # ---------------------------------------------------------------------
    # LÓGICA DE EXTRACCIÓN DE DATOS REAL
    # ---------------------------------------------------------------------
    try:
        print("INFO: Intentando conectar a la base de datos MySQL...")
        
        # 1. ESTABLECER CONEXIÓN 
        conn = mysql.connector.connect(
            host="localhost",        
            user="root",    
            password="arquisoft1", 
            database="dw_banco_de_datos" 
        )
        
        # 2. CONSULTA SQL (Tu consulta completa)
        SQL_QUERY = """
        SELECT 
            -- IDs
            v.Id_variante, 
            m.Id_muestra, 
            t.Id_tipo_tumor, 
            ug.Id_ubicacion_genomica,

            -- Datos de la muestra
            m.Tipo_de_muestra, 
            m.Fecha_obtencion, 
            m.Lugar_obtencion,

            -- Datos del tumor
            t.Nombre AS Tipo_tumor, 
            t.Descripcion AS Detalle_tumor, 
            t.Codigo_cie10,

            -- Datos de la variante
            v.Variante, 
            v.Tipo_modificacion, 
            v.Descripcion AS Detalle_variante,

            -- Datos de ubicación y clínica
            ug.Cromosoma, 
            ug.Brazo, 
            ug.Region,
            ug.Banda, 
            ug.Sub_banda,
            uxv.Significancia_clinica, 
            uxv.Frec_alelica

        FROM Variantes v
        LEFT JOIN Tipo_tumor_variante tv 
            ON v.Id_variante = tv.Id_variante
        LEFT JOIN Tipo_tumor t 
            ON tv.Id_tipo_tumor = t.Id_tipo_tumor
        LEFT JOIN Variante_x_informe vi 
            ON v.Id_variante = vi.Id_variante
        
        LEFT JOIN Informe_resultado i 
            ON vi.Id_informe = i.Id_informe
        LEFT JOIN Muestra_genomica m 
            ON i.Id_muestra = m.Id_muestra
        LEFT JOIN Ubicacion_x_variante uxv 
            ON v.Id_variante = uxv.Id_variante
        LEFT JOIN Ubicacion_genomica ug 
            ON uxv.Id_ubicacion_genomica = ug.Id_ubicacion_genomica
        ORDER BY m.Id_muestra, v.Id_variante;
        """
        
        # 3. EJECUTAR CONSULTA Y OBTENER DATOS EN UN DATAFRAME
        df = pd.read_sql(SQL_QUERY, conn)
        
        # 4. CERRAR CONEXIÓN
        conn.close()
        print(f"INFO: Conexión cerrada. Se extrajeron {len(df)} filas.")
        
        # ---------------------------------------------------------------------
        # EXPORTAR EL DATAFRAME A CSV
        # ---------------------------------------------------------------------
        df.to_csv(OUTPUT_CSV, index=False)
        print(f"✅ Éxito: Datos exportados a {OUTPUT_CSV}")

    except mysql.connector.Error as err:
        print(f"❌ ERROR MySQL: Algo salió mal al conectar o consultar: {err}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error inesperado durante la exportación: {e}")
        sys.exit(1)

# =========================================================================
# 3. EJECUCIÓN DEL SCRIPT (¡CORRECCIÓN APLICADA AQUÍ!)
# =========================================================================
if __name__ == "__main__":
    fetch_data_and_export()