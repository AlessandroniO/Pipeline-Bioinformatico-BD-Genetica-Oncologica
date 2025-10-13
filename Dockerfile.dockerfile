FROM ubuntu:22.04

# 1️⃣ Dependencias base y Miniconda
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y wget bzip2 ca-certificates git && \
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# 2️⃣ Configurar PATH y activar Conda base
ENV PATH="/opt/conda/bin:$PATH"
SHELL ["conda", "run", "-n", "base", "/bin/bash", "-c"]

# 3️⃣ Aceptar TOS
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# 4️⃣ Instalar mamba y limpiar
RUN conda install -n base -c conda-forge -y mamba && mamba clean --all -y

# 5️⃣ Instalar Snakemake y dependencias, incluyendo MySQL Connector
RUN mamba install -c bioconda -c conda-forge -y \
        snakemake \
        pandas \
        openpyxl \
        pyyaml \
        tqdm && \
    pip install mysql-connector-python && \
    mamba clean --all -y

# 6️⃣ Ajustes finales
WORKDIR /data
ENV CONDA_AUTO_UPDATE_CONDA=false
ENV PYTHONUNBUFFERED=1
ENV CONDA_SSL_VERIFY=false

# Activar base por defecto al iniciar bash
CMD ["conda", "run", "-n", "base", "bash"]
