import sys
import os
import pandas as pd
from sqlalchemy import create_engine
from dotenv import load_dotenv

# =========================================================================
# CONFIGURACIÓN Y CONEXIÓN
# =========================================================================

load_dotenv(os.path.join(os.path.dirname(__file__), '..', '..', '.env'))

DB_USER = os.getenv("MYSQL_USER")
DB_PASSWORD = os.getenv("MYSQL_ROOT_PASSWORD")
DB_HOST = os.getenv("MYSQL_HOST")
DB_PORT = os.getenv("MYSQL_PORT")
DB_NAME = os.getenv("MYSQL_DATABASE")

# =========================================================================
# FUNCIÓN AUXILIAR
# =========================================================================

def get_output_path():
    """Obtiene la ruta de salida desde argumentos o default."""
    if len(sys.argv) > 1:
        return sys.argv[1]
    return "data/processed/variantes_extraidas.csv"

# =========================================================================
# FUNCIÓN PRINCIPAL
# =========================================================================

def fetch_data_and_export():
    OUTPUT_CSV = get_output_path()
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)

    DB_URI = f"mysql+mysqlconnector://{DB_USER}:{DB_PASSWORD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"

    SQL_QUERY = """
            SELECT
                v.Id_variante AS id,
                v.Variante AS variante,
                v.Tipo_modificacion AS tipo_modificacion,
                v.Descripcion AS descripcion,
                v.Variante_id_clinvar as Id_clinvar,
                ug.Cromosoma AS chrom,
                ug.Brazo AS brazo,
                uxv.Frec_alelica AS frec_alelica,
                d.Barrio AS barrio,
                d.Ciudad AS ciudad,
                d.Provincia AS provincia,
                m.Tipo_de_muestra AS tipo_de_muestra
            FROM Variantes v
            LEFT JOIN Ubicacion_x_variante uxv ON v.Id_variante = uxv.Id_variante
            LEFT JOIN Ubicacion_genomica ug ON uxv.Id_ubicacion_genomica = ug.Id_ubicacion_genomica
            LEFT JOIN Variante_x_informe vxi ON v.Id_variante = vxi.Id_variante
            LEFT JOIN Informe_resultado ir ON vxi.Id_informe = ir.Id_informe
            LEFT JOIN Muestra_genomica m ON ir.Id_muestra = m.Id_muestra
            LEFT JOIN Paciente p ON m.Id_paciente = p.Id_paciente
            LEFT JOIN Direccion d ON p.Id_direccion = d.Id_direccion
            ORDER BY v.Id_variante;
            """

    try:
        engine = create_engine(DB_URI)
        df = pd.read_sql(SQL_QUERY, engine, dtype=str)


        # Limpieza: convertir a string y reemplazar valores vacíos por NA
        for col in df.columns:
            df[col] = df[col].astype(str).replace(['None', 'nan', 'NULL', 'NaN', 'N/A', 'NA'], pd.NA)

        print(f"INFO: Datos extraídos: {len(df)} filas.")

        # Agregar campo de origen
        df["source"] = "BD"

        # Exportar CSV
        df.to_csv(OUTPUT_CSV, index=False)
        print(f"✅ Éxito: Datos exportados a {OUTPUT_CSV} ({len(df)} filas)")

    except Exception as e:
        print(f"❌ ERROR durante la extracción: {e.__class__.__name__}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    fetch_data_and_export()
