# Snakemake Workflow para el TFG (Separación Somática/Germinal y QC)
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

# Archivo de salida de QC de Metadatos
METADATA_QC_REPORT = "data/results/qc/metadata_report.html"


# ----------------
# REGLA FINAL (TARGET)
# ----------------
# El objetivo final es generar los VCFs normalizados Y el reporte de QC.
rule all:
    input:
        FINAL_SOMATIC_VCF,
        FINAL_GERMLINE_VCF,
        METADATA_QC_REPORT


# ----------------------------------------------------
# ETAPA 0: ENRIQUECIMIENTO DE LA BASE DE DATOS (ETL)
# ----------------------------------------------------
rule db_enrichment:
    input: "GRCh38_latest_clinvar.vcf.gz"
    output: "data/intermediate/db_enriched.flag"
    conda: "envs/db_enrichment.yml"
    shell:
        "python pipeline/modules/etl_clinvar_variantes.py && touch {output}"


# ----------------------------------------------------
# ETAPA 1: EXTRACCIÓN DE DATOS DE LA BASE DE DATOS
# ----------------------------------------------------
rule db_extraction:
    input:
        # Dependencia agregada: debe esperar a que la BD esté enriquecida
        enriched_flag="data/intermediate/db_enriched.flag" 
    output:
        "data/intermediate/variantes_extraidas.csv"
    conda:
        "envs/db_extractor.yml"
    shell:
        "python pipeline/modules/db_extraction.py {output}"


# ------------------------------------------------------------------
# ETAPA 2A: PROCESAMIENTO Y DIVISIÓN DE DATOS (SOMÁTICO/GERMINAL)
# ------------------------------------------------------------------
rule data_processing_split:
    input:
        "data/intermediate/variantes_extraidas.csv"
    output:
        somatic_csv="data/processed/variantes_somaticas_limpias.csv",
        germline_csv="data/processed/variantes_germinales_limpias.csv"
    conda:
        "envs/processor.yml"
    shell:
        "python pipeline/modules/data_processing.py {input} {output.somatic_csv} {output.germline_csv}"


# --------------------------------------------------------------------
# NUEVA REGLA DE QC DE METADATOS (SE EJECUTA EN PARALELO)
# --------------------------------------------------------------------
rule metadata_qc_report:
    input:
        "data/intermediate/variantes_extraidas.csv"
    output:
        html="data/results/qc/metadata_report.html",
        errors_csv="data/results/qc/metadata_qc_errors.csv"
    conda:
        "envs/processor.yml" # Reutilizamos el entorno de Python/Pandas
    shell:
        "python pipeline/modules/metadata_qc.py {input} {output.html} {output.errors_csv}"


# --------------------------------------------------------------------
# ETAPA 2B: CONVERSIÓN A VCF - RUTA SOMÁTICA
# --------------------------------------------------------------------
rule vcf_conversion_somatic:
    input:
        "data/processed/variantes_somaticas_limpias.csv"
    output:
        "data/intermediate/somatic_raw.vcf"
    conda:
        "envs/processor.yml"
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
        "envs/processor.yml"
    shell:
        "python pipeline/modules/vcf_conversion.py {input} {output}"


# ----------------------------------------------------
# ETAPA 3A: NORMALIZACIÓN y COMPRESIÓN - SOMÁTICA
# ----------------------------------------------------
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
