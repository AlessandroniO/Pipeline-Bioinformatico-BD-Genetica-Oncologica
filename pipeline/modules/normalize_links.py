import pandas as pd
from sqlalchemy import create_engine
import os
import sys
from dotenv import load_dotenv

load_dotenv(os.path.join(os.path.dirname(__file__), '..', '..', '.env'))

DB_USER = os.getenv("MYSQL_USER")
DB_PASSWORD = os.getenv("MYSQL_ROOT_PASSWORD")
DB_HOST = os.getenv("MYSQL_HOST_LOCAL", "127.0.0.1") 
DB_PORT = os.getenv("MYSQL_PORT", "3307") 
DB_NAME = "dw_banco_de_datos" 
DB_URI = f"mysql+mysqlconnector://{DB_USER}:{DB_PASSWORD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"

# ESTO ES CRUCIAL:
INFO_VARIANTE_CSV = 'data/Info_variantes.csv' # Asegúrate de la ruta correcta
# Columna en info_variantes.csv que se usa para enlazar con la tabla Informe_resultado
# Debe ser un ID que exista en ambas fuentes, EJEMPLO: DNI del paciente.
CLAVE_DE_ENLACE = 'DNI_PACIENTE_O_ID_INFORME' # <<<<< DEBES AJUSTAR AQUÍ

def normalize_links(csv_path):
    try:
        engine = create_engine(DB_URI)
        df_variantes_info = pd.read_csv(csv_path)

        with engine.begin() as conn:
            
            # 1. Obtener IDs de tablas ya llenadas
            # Variante: (los 9781 que ya insertaste)
            df_variantes_ids = pd.read_sql("SELECT Id_variante, Variante FROM Variantes", conn)
            
            # Informe_resultado: Necesitamos el ID para el enlace
            # Para hacer el merge, necesitamos una columna común con info_variantes.csv
            df_informes_ids = pd.read_sql("SELECT ir.Id_informe, p.DNI FROM Informe_resultado ir JOIN Muestra_genomica m ON ir.Id_muestra = m.Id_muestra JOIN Paciente p ON m.Id_paciente = p.Id_paciente", conn)
            
            # 2. CREAR EL ENLACE CRUCIAL: VARIANTE_X_INFORME
            
            # 2a. Renombrar la clave de enlace del CSV de variantes
            df_enlaces = df_variantes_info.rename(columns={
                CLAVE_DE_ENLACE: 'DNI',  # Asumiendo que la clave de enlace es el DNI
                # Asumiendo que la variante está en una columna llamada 'Variante' o similar en info_variantes.csv
                'ID_DE_LA_VARIANTE_EN_EL_CSV': 'Variante' # <<<<< AJUSTAR AQUÍ si el nombre de la variante es diferente
            })

            # 2b. Merge con IDs de Variante (Para obtener Id_variante)
            df_enlaces = df_enlaces.merge(df_variantes_ids, on='Variante', how='inner')
            
            # 2c. Merge con IDs de Informe (Para obtener Id_informe)
            df_enlaces = df_enlaces.merge(df_informes_ids, on='DNI', how='inner') 
            
            # 2d. Insertar el enlace
            df_enlaces_final = df_enlaces[['Id_variante', 'Id_informe']].drop_duplicates()
            
            # Cargar Variante_x_informe
            df_enlaces_final.to_sql(
                'Variante_x_informe', conn, if_exists='append', index=False
            )
            print(f"INFO: {len(df_enlaces_final)} enlaces Variante-Informe creados (rellenando el eslabón roto).")
            
            # 3. Cargar Ubicacion_x_variante (Enriquecimiento Genómico)
            # Asumiendo que info_variantes.csv también tiene las columnas de enriquecimiento
            
            df_enrichment = df_variantes_info.rename(columns={
                'FRECUENCIA_ALELICA_CSV': 'Frecuencia_alelica', # <<<<< AJUSTAR
                'ID_UBICACION_GENOMICA_CSV': 'Id_ubicacion_genomica' # <<<<< AJUSTAR
            })
            
            # Este es un INSERT o UPDATE más complejo que omitiremos aquí, pero es el paso para llenar crom/brazo
            # df_enrichment[['Frecuencia_alelica', 'Id_variante', 'Id_ubicacion_genomica']].to_sql(
            #     'Ubicacion_x_variante', conn, if_exists='append', index=False
            # )
            # print("INFO: Datos de enriquecimiento de frecuencia/ubicación cargados.")

        print("✅ Éxito: Normalización de enlaces completada.")
        
    except Exception as e:
        print(f"❌ ERROR durante la normalización de enlaces: {e.__class__.__name__}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    normalize_links(INFO_VARIANTE_CSV)