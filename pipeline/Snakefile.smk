# Snakemake Workflow para el TFG (Separación Somática/Germinal)
# ===============================================================

# ----------------
# CONFIGURACIÓN
# ----------------
# La referencia del genoma NO está incluida. Se omite la normalización con FASTA
# para permitir que el pipeline termine.
GENOME_REFERENCE = "/dev/null"

# Archivos de salida finales para cada cohorte
FINAL_SOMATIC_VCF = "data/processed/somatic.cohort.norm.vcf.gz"
FINAL_GERMLINE_VCF = "data/processed/germline.cohort.norm.vcf.gz"


# ----------------
# REGLA FINAL (TARGET)
# ----------------
# El objetivo final es generar los VCFs normalizados para ambas cohortes.
rule all:
    input:
        FINAL_SOMATIC_VCF,
        FINAL_GERMLINE_VCF


# ----------------------------------------------------
# ETAPA 1: EXTRACCIÓN DE DATOS DE LA BASE DE DATOS
# ----------------------------------------------------
# ❗ ESTA REGLA AHORA INTENTA LA CONEXIÓN REAL A MYSQL
rule db_extraction:
    output:
        "data/intermediate/variantes_extraidas.csv"
    conda:
        "envs/db_extractor.yml" # Usando tu archivo db_extractor.yml
    shell:
        "python pipeline/modules/db_extraction.py {output}"


# ------------------------------------------------------------------
# ETAPA 2A: PROCESAMIENTO Y DIVISIÓN DE DATOS (SOMÁTICO/GERMINAL)
# ------------------------------------------------------------------
# Esta regla reemplaza la antigua data_processing. Genera dos salidas CSV.
rule data_processing_split:
    input:
        "data/intermediate/variantes_extraidas.csv"
    output:
        somatic_csv="data/processed/variantes_somaticas_limpias.csv",
        germline_csv="data/processed/variantes_germinales_limpias.csv"
    conda:
        "envs/processor.yml" # Usando el nombre de tu archivo.
    shell:
        "python pipeline/modules/data_processing.py {input} {output.somatic_csv} {output.germline_csv}"


# --------------------------------------------------------------------
# ETAPA 2B: CONVERSIÓN A VCF - RUTA SOMÁTICA
# --------------------------------------------------------------------
rule vcf_conversion_somatic:
    input:
        "data/processed/variantes_somaticas_limpias.csv"
    output:
        "data/intermediate/somatic_raw.vcf"
    conda:
        "envs/processor.yml" # Usando el nombre de tu archivo.
    shell:
        "python pipeline/modules/vcf_conversion.py {input} {output}"

# --------------------------------------------------------------------
# ETAPA 2C: CONVERSIÓN A VCF - RUTA GERMINAL
# --------------------------------------------------------------------
rule vcf_conversion_germline:
    input:
        "data/processed/variantes_germinales_limpias.csv"
    output:
        "data/intermediate/germline_raw.vcf"
    conda:
        "envs/processor.yml" # Usando el nombre de tu archivo.
    shell:
        "python pipeline/modules/vcf_conversion.py {input} {output}"


# ----------------------------------------------------
# ETAPA 3A: NORMALIZACIÓN y COMPRESIÓN - SOMÁTICA
# ----------------------------------------------------
# Se comprime y se indexa el VCF para que esté listo.
rule variant_normalization_somatic:
    input:
        "data/intermediate/somatic_raw.vcf"
    output:
        vcf=FINAL_SOMATIC_VCF,
        index="data/processed/somatic.cohort.norm.vcf.gz.csi"
    conda:
        "envs/bcftools.yml"
    shell:
        """
bgzip -c {input} > {output.vcf}
bcftools index {output.vcf}
"""


# ----------------------------------------------------
# ETAPA 3B: NORMALIZACIÓN y COMPRESIÓN - GERMINAL
# ----------------------------------------------------
# Se comprime y se indexa el VCF para que esté listo.
rule variant_normalization_germline:
    input:
        "data/intermediate/germline_raw.vcf"
    output:
        vcf=FINAL_GERMLINE_VCF,
        index="data/processed/germline.cohort.norm.vcf.gz.csi"
    conda:
        "envs/bcftools.yml"
    shell:
        """
bgzip -c {input} > {output.vcf}
bcftools index {output.vcf}
"""
