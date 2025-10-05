# =======================================================
# Archivo: pipeline/Snakefile
# =======================================================

rule all:
    input:
        "data/processed/variantes_extraidas.csv" 

rule db_extraction:
    # 1. Archivo de salida que esta regla producirá
    output:
        "data/processed/variantes_extraidas.csv"
    
    # 2. Archivo Conda (confirmamos el nombre correcto)
    conda:
        "envs/db_extractor.yml" 
        
    # 3. El comando shell con la ruta corregida a 'modules/'
    shell:
        "python pipeline/modules/export_variants2.py {output}" 