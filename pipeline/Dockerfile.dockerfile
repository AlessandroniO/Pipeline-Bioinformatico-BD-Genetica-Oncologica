FROM ubuntu:22.04

# =======================
# Configuración de Variables de Entorno
# =======================
# Evita prompts interactivos durante la instalación de paquetes.
ENV DEBIAN_FRONTEND=noninteractive
# Desactiva el buffering de Python (mejor para logs en contenedores).
ENV PYTHONUNBUFFERED=1
# Desactiva las actualizaciones automáticas de Conda.
ENV CONDA_AUTO_UPDATE_CONDA=false
# Precaución: Desactivar verificación SSL solo si es estrictamente necesario,
# por lo general, se recomienda mantenerlo habilitado.
# ENV CONDA_SSL_VERIFY=false

# =======================
# 1️⃣ Instalación de base, Miniconda y Mamba
# =======================
# Combinamos la instalación base, la descarga de Miniconda y la instalación de Mamba
# en una sola capa para optimizar el tamaño final de la imagen.
RUN apt-get update && \
    apt-get install -y wget bzip2 ca-certificates git && \
    \
    # Descargar e instalar Miniconda
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    \
    # Añadir Conda/Mamba al PATH temporal para el resto del comando RUN
    export PATH="/opt/conda/bin:${PATH}" && \
    \
    # 🚨 CORRECCIÓN: Aceptar los ToS antes de instalar mamba 🚨
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    \
    # Instalar Mamba desde conda-forge (más rápido para resolver dependencias)
    conda install -n base -c conda-forge -y mamba && \
    \
    # Limpieza de paquetes y Conda/Mamba para reducir el tamaño final
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    mamba clean --all -y

# Establecer el PATH final para todas las capas siguientes
ENV PATH="/opt/conda/bin:${PATH}"

# =======================
# 2️⃣ Instalación de dependencias del proyecto con Mamba
# =======================
# Instalamos Snakemake y las librerías de Python en una capa separada.
# Ahora incluye SQLAlchemy, necesario para las conexiones a bases de datos.
# -------------------------------------------------------------
# Paso 3: Instalar Snakemake y herramientas clave
# -------------------------------------------------------------
RUN mamba install -c bioconda -c conda-forge -y \
    snakemake && \
    mamba clean -y --all

# -------------------------------------------------------------
# Paso 4: Instalar librerías de Python (más rápido, sin bioconda)
# -------------------------------------------------------------
RUN mamba install -c conda-forge -y \
    cbioportal-client \
    pandas \
    sqlalchemy \
    biopython \
    psycopg2 \
    hgvs && \
    mamba clean -y --all

# =======================
# 3️⃣ Ajustes finales
# =======================
WORKDIR /data

