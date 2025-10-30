# Snakemake Workflow para el TFG (Separación Somática/Germinal y QC)
# ===============================================================

# ----------------
# CONFIGURACIÓN
# ----------------
GENOME_REFERENCE = "/dev/null"

# Archivos de salida finales para cada cohorte
FINAL_SOMATIC_VCF = "data/processed/somatic.cohort.norm.vcf.gz"
FINAL_GERMLINE_VCF = "data/processed/germline.cohort.norm.vcf.gz"

# Archivo de salida de QC de Metadatos
METADATA_QC_REPORT = "data/results/qc/metadata_report.html"

# Archivo de la lista de tumores CIE-10
CIE10_TUMORS = "Lista_tumores_CIE-10.csv"

# NUEVO: Salida de la consulta a cBioPortal
CBIOPORTAL_VARIANTS_CSV = "data/intermediate/cbioportal_frequent_variants.csv"


# ----------------
# REGLA FINAL (TARGET)
# ----------------
rule all:
    input:
        FINAL_SOMATIC_VCF,
        FINAL_GERMLINE_VCF,
        METADATA_QC_REPORT,
        CBIOPORTAL_VARIANTS_CSV # <-- NUEVA DEPENDENCIA DE SALIDA


# --------------------------------------------------------
# ETAPA 0.5: CONSULTA DE VARIANTES ONCOLÓGICAS (cBioPortal)
# --------------------------------------------------------
rule fetch_cbioportal_data:
    output:
        "data/intermediate/cbioportal_frequent_variants.csv"
    conda: "envs/fetch_cbioportal.yaml"

    shell:
        "python pipeline/modules/fetch_variants.py {output}"



# ----------------------------------------------------
# ETAPA 1: EXTRACCIÓN DE DATOS DE LA BASE DE DATOS
# ----------------------------------------------------
rule db_extraction:
    output:
        "data/processed/variantes_extraidas.csv"
    conda:
        "envs/db_extractor.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        python pipeline/modules/db_extraction.py {output}
        """


# Variables globales / Archivos estáticos
# o que la variable global CIE10_TUMORS apunta a la ruta donde se encuentra el archivo que subiste.
CIE10_TUMORS = "Lista_tumores_CIE-10.csv" 
# NOTA: Asegúrate de que esta ruta sea correcta, o usa el nombre del archivo si está en el directorio raíz:
# CIE10_TUMORS = "Lista_tumores_CIE-10.csv" 


# Variables globales / Archivos estáticos
# NOTA: Ajusta la ruta a tu archivo CIE-10 (asumo que es un archivo de configuración estático).
CIE10_TUMORS = "Lista_tumores_CIE-10.csv" 

# ------------------------------------------------------------------
# ETAPA 2: PROCESAMIENTO Y DIVISIÓN DE DATOS (SOLO DB y CIE-10)
# ------------------------------------------------------------------
rule data_processing_split:
    input:
        db_data="data/processed/variantes_extraidas.csv", # Output de db_extraction
        cie10_map=CIE10_TUMORS 
    output:
        somatic_csv="data/processed/variantes_somaticas_limpias.csv",
        germline_csv="data/processed/variantes_germinales_limpias.csv",
        qc_summary="data/results/qc/qc_summary.csv" 
    conda:
        "envs/processor.yml"
    shell:
        # Pasamos solo los 4 argumentos necesarios: db_csv, cie10_file, somatic_out, germline_out
        "python pipeline/modules/data_processing.py {input.db_data} "
        "{input.cie10_map} {output.somatic_csv} {output.germline_csv}"

# --------------------------------------------------------------------
# REGLA DE QC DE METADATOS (SE EJECUTA EN PARALELO)
# --------------------------------------------------------------------
rule metadata_qc_report:
    input:
        "data/intermediate/variantes_extraidas.csv"
    output:
        html="data/results/qc/metadata_report.html",
        errors_csv="data/results/qc/metadata_qc_errors.csv"
    conda:
        "envs/processor.yml"
    shell:
        "python pipeline/modules/metadata_qc.py {input} {output.html} {output.errors_csv}"


# ======================================================================
# Variables de Ruta Corregidas (Tus archivos VCF y de intermedios)
# ======================================================================
# La ruta real de tu VCF de referencia
CLINVAR_VCF_FILE = "data/raw/GRCh38_latest_clinvar.vcf.gz.vcf" 

# Archivo intermedio que guardará el VCF parseado a TSV
PARSED_VCF_TSV = "data/intermediate/clinvar_grch38_parsed.tsv" 

rule convert_somatic_to_genomic:
    input:
        csv_in="data/processed/variantes_somaticas_limpias.csv"
    output:
        csv_out="data/processed/variantes_somaticas_genomic.csv"
    conda:
        "envs/hgvs.yml"
    shell:
        """
        python pipeline/modules/convert_to_genomics.py {input.csv_in} {output.csv_out}
        """

rule convert_germline_to_genomic:
    input:
        csv_in="data/processed/variantes_germinales_limpias.csv"
    output:
        csv_out="data/processed/variantes_germinales_genomic.csv"
    conda:
        "envs/hgvs.yml"
    shell:
        """
        python pipeline/modules/convert_to_genomics.py {input.csv_in} {output.csv_out}
        """


# ======================================================================
# Regla 1: Pre-procesar el VCF (VCF a TSV) - NUEVA REGLA
# ======================================================================
rule vcf_to_tsv:
    input:
        vcf_in=CLINVAR_VCF_FILE 
    output:
        tsv_out=PARSED_VCF_TSV 
    conda:
        "envs/processor.yml" 
    shell:
        """
        python pipeline/modules/vcf_parser.py {input.vcf_in} {output.tsv_out}
        """

# ======================================================================
# Regla 2.A: ANOTACIÓN POR JOIN - SOMÁTICA (SUSTITUYE annotate_api_somatic)
# ======================================================================
# Usamos un nuevo nombre interno (annotate_join_somatic) para la regla,
# pero mantenemos el nombre de salida (variantes_somaticas_API_anotadas.csv).
rule annotate_join_somatic:
    input:
        csv_in="data/processed/variantes_somaticas_genomic.csv",
        ref_tsv=PARSED_VCF_TSV # Depende de la regla vcf_to_tsv
    output:
        # CRÍTICO: Mantenemos el nombre de archivo de la antigua API.
        csv_out="data/intermediate/variantes_somaticas_API_anotadas.csv" 
    conda:
        "envs/processor.yml" 
    shell:
        """
        python pipeline/modules/join_with_vcf.py {input.csv_in} {input.ref_tsv} {output.csv_out}
        """

# ======================================================================
# Regla 2.B: ANOTACIÓN POR JOIN - GERMINAL (SUSTITUYE annotate_api_germline)
# ======================================================================
# Usamos un nuevo nombre interno (annotate_join_germline) para la regla,
# pero mantenemos el nombre de salida (variantes_germinales_API_anotadas.csv).
rule annotate_join_germline:
    input:
        csv_in="data/processed/variantes_germinales_genomic.csv",
        ref_tsv=PARSED_VCF_TSV 
    output:
        # CRÍTICO: Mantenemos el nombre de archivo de la antigua API.
        csv_out="data/intermediate/variantes_germinales_API_anotadas.csv"
    conda:
        "envs/processor.yml"
    shell:
        """
        python pipeline/modules/join_with_vcf.py {input.csv_in} {input.ref_tsv} {output.csv_out}
        """


