import sys
import pandas as pd

def prepare_vep_input(input_csv, output_txt):
    """
    Lee el CSV de variantes extraídas, extrae las variantes únicas de la
    columna 'variante' (notación HGVS) y las guarda en un archivo de texto
    simple, una variante por línea, para usar como entrada de VEP.
    """
    try:
        # Cargar el CSV
        df = pd.read_csv(input_csv)
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo de entrada: {input_csv}", file=sys.stderr)
        sys.exit(1)

    # 1. Extraer la notación HGVS de nucleótido
    # La columna 'variante' tiene el formato: NM_000321.3(RB1):c.2106+22C>G (p.Val...)
    # VEP solo necesita la parte de la variante (ej: NM_000321.3:c.2106+22C>G)
    
    # Creamos una columna temporal quitando la parte de proteína (p....)
    df['hgvs_clean'] = df['variante'].str.split(' ').str[0]
    
    # 2. Extraer las variantes únicas
    # Es crucial trabajar solo con las variantes únicas para no sobrecargar VEP.
    unique_variants = df['hgvs_clean'].unique()

    # 3. Guardar en el formato de entrada de VEP
    print(f"Escribiendo {len(unique_variants)} variantes únicas para VEP en {output_txt}")
    with open(output_txt, 'w') as f:
        # VEP espera una variante por línea
        f.write('\n'.join(unique_variants))

    print("Archivo de entrada de VEP preparado con éxito.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python prepare_vep_input.py <input_csv> <output_txt>")
        sys.exit(1)

    input_csv = sys.argv[1]
    output_txt = sys.argv[2]
    
    prepare_vep_input(input_csv, output_txt)