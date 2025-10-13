FROM ubuntu:22.04

# =======================
# 1️⃣ Instalar dependencias base y Miniconda
# =======================
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y wget bzip2 ca-certificates git && \
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV PATH="/opt/conda/bin:${PATH}"

# =======================
# 2️⃣ Aceptar los Términos de Servicio (TOS)
# =======================
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Paso 3️⃣ Instalar mamba (rápido)
RUN conda install -n base -c conda-forge -y mamba && \
    mamba clean --all -y

# Paso 4️⃣ Instalar snakemake (en paso aparte, menos propenso a fallos)
RUN mamba install -c bioconda -c conda-forge -y snakemake && \
    mamba clean --all -y

# =======================
# 5️⃣ Ajustes finales
# =======================
WORKDIR /data
ENV CONDA_AUTO_UPDATE_CONDA=false
ENV PYTHONUNBUFFERED=1

# Evita errores de conexión durante el build
ENV CONDA_SSL_VERIFY=false
