# 🧬 Pipeline Bioinformático para el Banco de Datos Genómico-Oncológico

Este repositorio contiene un **pipeline bioinformático reproducible** diseñado para la **extracción, clasificación, anotación y normalización de variantes genéticas** dentro del marco del **Banco de Datos Genómico-Oncológico (BDGO)**, desarrollado en el CIBIOIMP – Universidad Católica de Córdoba.

---

## 📖 Descripción general

El pipeline permite integrar información genómica proveniente de diferentes fuentes (bases de datos, archivos clínicos o secuenciación) para generar **cohortes somáticas y germinales** normalizadas y listas para análisis posteriores.  

Incluye pasos de:
- Extracción y limpieza de datos desde la base relacional MySQL.
- Clasificación de variantes según evidencia **COSMIC**.
- Anotación cruzada con **ClinVar**, **gnomAD** y bases de referencia.
- Normalización y escritura en formato **VCF** estandarizado.
- Generación de resúmenes, métricas y reportes de control de calidad.

Cada etapa se ejecuta automáticamente mediante **Snakemake**, garantizando reproducibilidad, modularidad y trazabilidad completa.

---

## ⚙️ Estructura del repositorio

Pipeline-Bioinformatico-BD-Genetica-Oncologica/
├── data/
│ ├── processed/ # Resultados intermedios (ignorado en Git)
│ ├── intermediate/ # Archivos de trabajo y temporales (ignorado)
│ └── results/ # Salidas finales (tablas, CSV, QC)
├── pipeline/
│ ├── Snakefile.smk # Definición principal del flujo Snakemake
│ ├── Dockerfile.dockerfile# Imagen base del pipeline
│ ├── envs/ # Entornos Conda para cada módulo
│ └── modules/ # Scripts Python de cada etapa
│ ├── db_extraction.py
│ ├── cosmic_classify.py
│ ├── annotate_with_gnomad.py
│ ├── merge_gnomad_back.py
│ ├── csv_to_vcf.py
│ ├── make_final_summaries.py
│ └── ...
├── .gitignore
└── README.md

## 🚀 Requisitos

- **Sistema operativo:** Ubuntu ≥ 22.04 LTS o Windows WSL2  
- **Snakemake** ≥ 8.0  
- **Conda / Mamba** ≥ 1.5  
- **Docker** (opcional, para ejecución contenida)  
- **Python** ≥ 3.11


## 🧩 Instalación

```bash
# Clonar el repositorio
git clone https://github.com/AlessandroniO/Pipeline-Bioinformatico-BD-Genetica-Oncologica.git
cd Pipeline-Bioinformatico-BD-Genetica-Oncologica

# Crear y activar entorno base de Snakemake
mamba create -n snakemake -c bioconda snakemake
conda activate snakemake
▶️ Ejecución del pipeline
Modo local (Conda)
bash
Copiar código
snakemake --cores 4 --use-conda
Modo Docker (recomendado)
bash
Copiar código
docker build -t bioinfo-pipeline:latest -f pipeline/Dockerfile.dockerfile .
docker run --rm -v $(pwd):/app -w /app bioinfo-pipeline:latest snakemake --cores 4
Cada regla del pipeline ejecuta un módulo Python independiente y deja su log correspondiente en logs/.

📊 Etapas principales
Etapa	Script / Módulo	Descripción
1️⃣ Extracción	db_extraction.py	Obtiene variantes de la base de datos relacional MySQL.
2️⃣ Clasificación	cosmic_classify.py	Clasifica variantes como somáticas o germinales según evidencia COSMIC.
3️⃣ Anotación	annotate_with_gnomad.py, merge_gnomad_back.py	Enlaza variantes con frecuencias poblacionales y anotaciones de ClinVar / gnomAD.
4️⃣ Normalización	csv_to_vcf.py	Convierte los resultados a formato VCF estándar.
5️⃣ Reportes	make_final_summaries.py	Genera tablas de resumen y reportes de control de calidad.

🧠 Filosofía de diseño
Reproducible: todas las dependencias se manejan con Conda y Docker.

Modular: cada etapa se puede ejecutar o reutilizar de forma independiente.

Transparente: logs, rutas y parámetros configurables desde config/ y Snakefile.smk.

Interoperable: salidas compatibles con VCFtools, bcftools, ClinVar, gnomAD, Power BI, etc.

⚖️ Licencia y atribución
Este proyecto se distribuye bajo licencia MIT.
Desarrollado por Antonio María Garzón como parte del Trabajo Final de Grado de la Licenciatura en Bioinformática – Universidad Católica de Córdoba (UCC), en el marco del Centro de Investigación en Bioinformática y Medicina de Precisión (CIBIOIMP).

🧩 Contacto
📧 2115318@ucc.edu.ar
🏛️ CIBIOIMP – Facultad de Ciencias Químicas – Universidad Católica de Córdoba
🌐 https://github.com/AlessandroniO
