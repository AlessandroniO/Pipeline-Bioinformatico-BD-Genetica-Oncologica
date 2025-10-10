import sys
import os
import pandas as pd

# =========================================================================
# 1. FUNCIÓN PRINCIPAL DE QC DE METADATOS
# =========================================================================
def run_metadata_qc(input_csv_path, output_html_path, output_errors_path):
    """
    Realiza el Control de Calidad (QC) de las columnas de metadatos clave 
    (tipo_tumor, sexo, edad, barrio) y genera un reporte HTML y un CSV de errores.
    """
    
    # Crear el directorio de resultados
    os.makedirs(os.path.dirname(output_html_path), exist_ok=True)
    
    try:
        # Leer los datos.
        print(f"INFO: Leyendo datos extraídos para QC de metadatos desde: {input_csv_path}")
        df = pd.read_csv(input_csv_path)
    except Exception as e:
        print(f"❌ ERROR: No se pudo leer el CSV de entrada para QC: {e}")
        sys.exit(1)

    # ❗ CORRECCIÓN CRUCIAL: Normalizar todas las cabeceras a minúsculas
    df.columns = df.columns.str.lower()

    # Definir las columnas clave para el QC (ahora en minúsculas)
    QC_COLS = ['tipo_tumor', 'sexo', 'edad', 'barrio']
    
    # ---------------------------------------------------------------------
    # 2. VALIDACIÓN DE REGLAS
    # ---------------------------------------------------------------------
    qc_results = []
    error_rows = []

    # A. Detección de Valores Nulos
    for col in QC_COLS:
        # Verificar que la columna existe después de la normalización
        if col not in df.columns:
            qc_results.append({
                'Campo': col,
                'Tipo_Error': 'Columna Faltante',
                'Total_Errores': len(df),
                'Detalle': f'La columna clave "{col}" no se encontró en el CSV.'
            })
            continue # Saltar a la siguiente columna si falta una
            
        nulos = df[col].isnull().sum()
        qc_results.append({
            'Campo': col,
            'Tipo_Error': 'Valores Nulos (NaN)',
            'Total_Errores': nulos,
            'Detalle': f'Se encontraron {nulos} valores nulos/vacíos.'
        })
        if nulos > 0:
            error_rows.append(df[df[col].isnull()].assign(Error_QC=f"Nulo en {col}"))

    # B. Validación de Sexo (solo 'M', 'F', 'Otro' o nulo)
    sexo_col = 'sexo'
    if sexo_col in df.columns:
        valid_sexo = ['m', 'f']
        # Buscamos valores que no estén en la lista de válidos Y que no sean nulos
        errores_sexo = df[~df[sexo_col].astype(str).str.lower().isin(valid_sexo) & df[sexo_col].notnull()]
        n_errores_sexo = len(errores_sexo)
        
        qc_results.append({
            'Campo': sexo_col,
            'Tipo_Error': 'Categoría Inválida',
            'Total_Errores': n_errores_sexo,
            'Detalle': f"Categorías fuera de ['M', 'F']. Ejemplos: {errores_sexo[sexo_col].unique()[:3].tolist()}"
        })
        if n_errores_sexo > 0:
            error_rows.append(errores_sexo.assign(Error_QC="Valor de Sexo Inválido"))

    # C. Validación de Edad (fuera de rango lógico: 0-100)
    edad_col = 'edad'
    if edad_col in df.columns:
        # Convertir a numérico, forzando errores a NaN
        edades_validas = pd.to_numeric(df[edad_col], errors='coerce')
        
        # Errores en rango (y asegurando que no son los nulos ya contados)
        errores_edad_rango = df[((edades_validas < 0) | (edades_validas > 100)) & (df[edad_col].notnull())]
        n_errores_edad = len(errores_edad_rango)
        
        qc_results.append({
            'Campo': edad_col,
            'Tipo_Error': 'Rango Inválido (Edad)',
            'Total_Errores': n_errores_edad,
            'Detalle': f'Edades fuera del rango lógico (0-100 años).'
        })
        if n_errores_edad > 0:
            error_rows.append(errores_edad_rango.assign(Error_QC="Edad fuera de 0-100"))

    # D. Detección de 'Barrios' raros (categorías con menos de 2 ocurrencias)
    barrio_col = 'barrio'
    if barrio_col in df.columns:
        conteo_barrios = df[barrio_col].value_counts()
        barrios_raros = conteo_barrios[conteo_barrios < 2].index
        
        n_barrios_raros = len(df[df[barrio_col].isin(barrios_raros)])
        
        qc_results.append({
            'Campo': barrio_col,
            'Tipo_Error': 'Categoría Rara/Única',
            'Total_Errores': n_barrios_raros,
            'Detalle': f'Barrios con menos de 2 ocurrencias. Total de {len(barrios_raros)} categorías únicas a revisar.'
        })

    # ---------------------------------------------------------------------
    # 3. GENERACIÓN DE REPORTES
    # ---------------------------------------------------------------------
    
    # Concatenar filas de error y guardar CSV para revisión manual
    if error_rows:
        final_errors_df = pd.concat(error_rows, ignore_index=True)
        final_errors_df.to_csv(output_errors_path, index=False)
        print(f"🚨 ADVERTENCIA: Se guardaron {len(final_errors_df)} filas con errores de QC en {output_errors_path}")
    else:
        # Crear un archivo de errores vacío si no hay problemas
        pd.DataFrame().to_csv(output_errors_path, index=False)
        print(f"✅ Éxito: No se detectaron errores críticos de QC. CSV de errores generado vacío en {output_errors_path}")


    # Generar Reporte HTML
    reporte_df = pd.DataFrame(qc_results)
    total_errores = reporte_df['Total_Errores'].sum()
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="es">
    <head>
        <meta charset="UTF-8">
        <title>Reporte de Control de Calidad (QC) de Metadatos</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; background-color: #f4f4f9; color: #333; }}
            .container {{ background-color: #fff; padding: 25px; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.1); }}
            h1 {{ color: #007bff; border-bottom: 2px solid #eee; padding-bottom: 10px; }}
            h2 {{ color: #ff5722; }}
            table {{ width: 100%; border-collapse: collapse; margin-top: 20px; }}
            th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
            th {{ background-color: #f0f0f0; }}
            .summary-box {{ padding: 15px; margin-bottom: 20px; border-radius: 4px; }}
            .success {{ background-color: #e6ffe6; border-left: 5px solid #4CAF50; }}
            .warning {{ background-color: #fff8e1; border-left: 5px solid #FFC107; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>Reporte de Control de Calidad de Metadatos</h1>
            <p>Generado por: <code>metadata_qc.py</code></p>
            <p>Datos de entrada: <code>{os.path.basename(input_csv_path)}</code></p>
            <div class="summary-box {'success' if total_errores == 0 else 'warning'}">
                <strong>Resumen:</strong> Se detectaron <strong>{total_errores}</strong> problemas de calidad en las columnas clave (tipo_tumor, sexo, edad, barrio).
            </div>
            
            <h2>Detalle de Errores por Campo</h2>
            {reporte_df.to_html(index=False)}

            <h2>Acción Requerida</h2>
            <p>Para la revisión manual y corrección, el archivo CSV completo con las filas problemáticas ha sido guardado en:</p>
            <code>{output_errors_path}</code>
        </div>
    </body>
    </html>
    """
    
    with open(output_html_path, 'w') as f:
        f.write(html_content)
        
    print(f"✅ Éxito: Reporte HTML de QC generado en {output_html_path}")

# =========================================================================
# 4. EJECUCIÓN
# =========================================================================
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python metadata_qc.py <input_csv> <output_html> <output_errors_csv>")
        sys.exit(1)
    
    input_csv = sys.argv[1]
    output_html = sys.argv[2]
    output_errors = sys.argv[3]
    
    run_metadata_qc(input_csv, output_html, output_errors)
