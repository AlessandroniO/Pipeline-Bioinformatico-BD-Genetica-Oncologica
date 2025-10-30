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

RAW_CLINICAL_CSV = 'data/raw/Datos_ficticios_pacientes_2012-2024.csv' 
# Ajusta los nombres de las columnas en este script a tus nombres reales del CSV clínico.

def normalize_clinical_data(csv_path):
    try:
        engine = create_engine(DB_URI)
        df = pd.read_csv(csv_path)

        with engine.begin() as conn:
            
            # --- OBTENER ID DE DIRECCION ---
            df_dir_ids = pd.read_sql("SELECT Id_direccion, Barrio, Ciudad, Provincia FROM Direccion", conn)
            
            # 1. NORMALIZAR E INSERTAR PACIENTES (AJUSTA LOS left_on)
            df_paciente = df.drop_duplicates(subset=['DNI_paciente']).rename(columns={
                'DNI_paciente': 'DNI', 'Nombre_paciente': 'Nombre', 
                'Apellido_paciente': 'Apellido', 'Fecha_nac_paciente': 'Fecha_nacimiento',
                'Sexo_paciente': 'Sexo', 'Obra_social_paciente': 'Obra_social', 
                'Nacionalidad_paciente': 'Nacionalidad', 'estado_vital': 'Estado_vital',
                'fecha_defuncion': 'Fecha_defuncion', 'causa_muerte': 'Causa_muerte',
                # Ajusta estas 3 columnas de tu CSV a los nombres reales
                'Direccion_paciente_barrio': 'Barrio_CSV', 
                'Direccion_paciente_ciudad': 'Ciudad_CSV', 
                'Direccion_paciente_provincia': 'Provincia_CSV',
            })
            
            df_paciente = df_paciente.merge(
                df_dir_ids,
                left_on=['Barrio_CSV', 'Ciudad_CSV', 'Provincia_CSV'],
                right_on=['Barrio', 'Ciudad', 'Provincia'],
                how='left'
            )
            df_paciente[['Nombre', 'Apellido', 'DNI', 'Fecha_nacimiento', 'Sexo', 'Obra_social', 
                         'Nacionalidad', 'Estado_vital', 'Fecha_defuncion', 'Causa_muerte', 'Id_direccion']].to_sql(
                'Paciente', conn, if_exists='append', index=False
            )
            print(f"INFO: {len(df_paciente)} pacientes insertados/actualizados.")
            
            # --- OBTENER IDs DE PACIENTE ---
            df_pacientes_ids = pd.read_sql("SELECT Id_paciente, DNI FROM Paciente", conn)
            
            # 2. NORMALIZAR E INSERTAR MUESTRA GENÓMICA
            df_muestra = df.drop_duplicates(subset=['DNI_paciente', 'fecha_obtencion_muestra']).rename(columns={
                'DNI_paciente': 'DNI', 'tipo_muestra': 'Tipo_de_muestra', 
                'fecha_obtencion_muestra': 'Fecha_obtencion', 'lugar_obtencion': 'Lugar_obtencion',
            })
            df_muestra = df_muestra.merge(df_pacientes_ids, on='DNI', how='left') 
            df_muestra[['Tipo_de_muestra', 'Fecha_obtencion', 'Lugar_obtencion', 'Id_paciente']].to_sql(
                'Muestra_genomica', conn, if_exists='append', index=False
            )
            print(f"INFO: {len(df_muestra)} muestras insertadas.")
            
            # --- OBTENER IDs DE MUESTRA ---
            df_muestras_ids = pd.read_sql("SELECT Id_muestra, Id_paciente FROM Muestra_genomica", conn)
            
            # 3. NORMALIZAR E INSERTAR INFORME RESULTADO
            df_informe = df.drop_duplicates(subset=['DNI_paciente', 'tipo_analisis']).rename(columns={
                'empresa_a_cargo': 'Empresa_encargada', 'tipo_analisis': 'Tipo_analisis',
            })
            df_informe = df_informe.merge(df_pacientes_ids, left_on='DNI_paciente', right_on='DNI', how='left')
            df_informe = df_informe.merge(df_muestras_ids, on='Id_paciente', how='left') # Enlazar con Id_muestra

            df_informe[['Empresa_encargada', 'Tipo_analisis', 'Id_muestra']].to_sql(
                'Informe_resultado', conn, if_exists='append', index=False
            )
            print(f"INFO: {len(df_informe)} informes insertados.")

        print("✅ Éxito: Normalización de datos clínicos completada.")
        
    except Exception as e:
        print(f"❌ ERROR durante la normalización clínica: {e.__class__.__name__}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    normalize_clinical_data(RAW_CLINICAL_CSV)