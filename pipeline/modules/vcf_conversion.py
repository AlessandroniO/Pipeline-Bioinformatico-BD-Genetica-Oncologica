import sys
import os
import pandas as pd

# =========================================================================
# 1. FUNCIÓN PRINCIPAL
# =========================================================================
def create_vcf(input_csv_path, output_vcf_path):
    """
    Lee un CSV de variantes limpias y genera un archivo VCF estándar (v4.3).
    
    Args:
        input_csv_path (str): Ruta al archivo CSV de entrada (somático o germinal).
        output_vcf_path (str): Ruta al archivo VCF de salida.
    """
    
    # 1. Crear el directorio de salida si no existe
    output_dir = os.path.dirname(output_vcf_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    try:
        # 2. Leer los datos y verificar la columna crucial (corregida a 'frec_alelica')
        print(f"INFO: Leyendo datos limpios de: {input_csv_path}")
        df = pd.read_csv(input_csv_path)
        
        # 3. VERIFICACIÓN DE COLUMNA CORREGIDA
        if 'frec_alelica' not in df.columns:
            print(f"❌ FATAL ERROR: No se puede acceder a la columna clave: 'frec_alelica'")
            print(f"🚨 Las columnas disponibles son: {list(df.columns)}")
            print("POSIBLE CAUSA: El script anterior (data_processing.py) no generó la columna correctamente.")
            sys.exit(1)

        # 4. Abrir el archivo VCF de salida
        with open(output_vcf_path, 'w') as vcf_file:
            
            # 5. Escribir el encabezado VCF (meta-información)
            vcf_file.write("##fileformat=VCFv4.3\n")
            vcf_file.write(f"##source=DW_Extraction_Pipeline\n")
            vcf_file.write(f'##INFO=<ID=DBID,Number=1,Type=String,Description="ID de la variante en la base de datos: Id_variante">\n')
            vcf_file.write(f'##INFO=<ID=TUMOR,Number=1,Type=String,Description="Tipo de tumor asociado">\n')
            vcf_file.write(f'##INFO=<ID=CS,Number=1,Type=String,Description="Significancia clínica">\n')
            vcf_file.write(f'##INFO=<ID=AF,Number=A,Type=Float,Description="Frecuencia Alélica (frec_alelica)">\n')
            
            # 6. Escribir la línea de cabecera de las columnas VCF
            vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

            # 7. Iterar sobre las filas del DataFrame y escribir el cuerpo del VCF
            for index, row in df.iterrows():
                
                # Campos obligatorios
                chrom = str(row['cromosoma']).replace('chr', '')
                # Usamos la variante para determinar REF y ALT.
                # Asumiendo que 'variante' tiene el formato [POS] [REF]>[ALT] o similar,
                # o usamos una lógica de conversión simple.
                # Dado que no tenemos POS exacto, usaremos Id_ubicacion_genomica como POS (temporalmente)
                # y simplificaremos REF/ALT hasta que se defina la columna VCF.
                
                # USANDO ID DE VARIANTE Y VARIACIÓN
                # Para simplificar la conversión, usaremos una posición y la variación.
                # Nota: Una conversión VCF real requiere columnas REF y ALT precisas.
                
                # Lógica simplificada basada en tu SQL:
                
                # Asumimos que la columna 'brazo' es lo suficientemente precisa para un POS de prueba.
                pos = str(row['id_ubicacion_genomica']) # Usamos el ID como POS temporal
                ref = "N" # Base de referencia genérica
                alt = str(row['variante']).split('>')[-1].strip() if '>' in str(row['variante']) else "V" # Base variante genérica

                # Si el cromosoma no tiene un formato válido, saltar la fila
                if not chrom or not pos or not alt:
                     print(f"ADVERTENCIA: Saltando variante {row['id_variante']} por datos incompletos o inválidos.")
                     continue
                
                # Campos INFO
                info_list = []
                info_list.append(f"DBID={row['id_variante']}")
                info_list.append(f"TUMOR={row['tipo_tumor']}")
                info_list.append(f"CS={row['significancia_clinica']}")
                
                # Corrección crucial: Usar 'frec_alelica' y formatearla
                af_val = row['frec_alelica'] if pd.notna(row['frec_alelica']) else 0.0
                info_list.append(f"AF={af_val:.4f}")

                INFO = ";".join(info_list)
                
                # Campos de formato (simples para este prototipo)
                FORMAT = "GT:AF"
                SAMPLE = f"1/1:{af_val:.4f}" # Genotipo 1/1 y AF
                
                # Línea VCF
                vcf_line = [
                    chrom,                  # CHROM
                    pos,                    # POS
                    row['id_variante'],     # ID
                    ref,                    # REF
                    alt,                    # ALT
                    ".",                    # QUAL (No definido)
                    "PASS",                 # FILTER
                    INFO,                   # INFO
                    FORMAT,                 # FORMAT
                    row['id_muestra']       # SAMPLE (Usamos el ID de muestra como nombre de la columna)
                ]
                
                vcf_file.write("\t".join(map(str, vcf_line)) + "\n")
        
        print(f"✅ Éxito: Conversión completada. VCF generado en {output_vcf_path}")

    except Exception as e:
        print(f"❌ ERROR CRÍTICO durante la conversión VCF: {e}")
        sys.exit(1)

# =========================================================================
# 2. EJECUCIÓN
# =========================================================================
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python vcf_conversion.py <input_csv> <output_vcf>")
        sys.exit(1)
    
    input_csv = sys.argv[1]
    output_vcf = sys.argv[2]
    
    create_vcf(input_csv, output_vcf)