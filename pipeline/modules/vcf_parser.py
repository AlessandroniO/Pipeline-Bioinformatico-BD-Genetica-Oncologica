import pandas as pd
import sys

def parse_vcf_to_tsv(input_vcf, output_tsv):
    """
    Lee el VCF y extrae las coordenadas VCF y el CAMPO INFO COMPLETO 
    para permitir la unión por subcadena HGVS en el siguiente paso.
    """
    print(f"Iniciando PARSING VCF para unión por SUBSTRING: {input_vcf}...")
    
    data_rows = []
    header = None
    variant_count = 0
    
    try:
        with open(input_vcf, 'r') as f:
            for line in f:
                if line.startswith('##'):
                    continue
                
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    continue
                
                if header is None or not line.strip(): 
                    continue

                fields = line.strip().split('\t')
                
                # Verificación de columnas mínimas (deberían ser 8 en un VCF estándar)
                if len(fields) < 8:
                    continue

                row_data = dict(zip(header, fields))
                
                # NO FILTRAMOS. GUARDAMOS TODA LA VARIANTE
                if row_data['ID'] and row_data['REF']:
                    data_rows.append({
                        '#CHROM': row_data['#CHROM'],
                        'POS': int(row_data['POS']),
                        'REF': row_data['REF'],
                        'ALT': row_data['ALT'],
                        # Guardamos el string INFO COMPLETO para buscar el HGVS después
                        'INFO_STRING': row_data['INFO']
                    })
                    variant_count += 1

    except Exception as e:
        print(f"ERROR CRÍTICO EN PARSER: Fallo al procesar el archivo VCF. {e}", file=sys.stderr)
        sys.exit(1)

    df_clean = pd.DataFrame(data_rows)
    df_clean.to_csv(output_tsv, sep='\t', index=False)
    
    print(f"PARSING COMPLETO. {variant_count} variantes guardadas en {output_tsv} con INFO_STRING.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python vcf_parser.py <input_vcf> <output_tsv>", file=sys.stderr)
        sys.exit(1)

    parse_vcf_to_tsv(sys.argv[1], sys.argv[2])