# alzheimer-gwas-curation

Curation, standardisation and overlap analysis of metadata and genotyping
platforms across forty Alzheimer's disease GWAS.

This repository contains the curation pipeline and analysis code developed during
a curricular internship at i3S (Instituto de Investigação e Inovação em Saúde), in
support of a methodological review of genome-wide association studies (GWAS) in
Alzheimer's disease. The work quantifies two forms of overlap between studies, at
the level of source cohorts and at the level of genetic markers (SNPs), and
characterises the heterogeneity of the genotyping platforms used.

## Overview

Genome-wide association studies (GWAS) of Alzheimer's disease have grown
considerably over the past two decades, yet larger samples have not translated
into a proportional gain in genetic signal. Two reasons stand out, the reuse of
the same source cohorts across studies, which undermines their independence, and
the heterogeneity of the genetic data, which come from different genotyping
platforms that are rarely harmonised. This project provides a quantitative basis
for both issues, turning into numbers what the literature had described mostly in
qualitative terms.

The work is organised around three strands:

- **Metadata curation.** A structured metadata table covering forty Alzheimer's
  disease GWAS was built from scratch, drawing together information on source
  cohorts, case and control counts, genomic assay type, diagnostic methods and
  reported cofactors. This table is the source on which every later analysis
  depends.
- **Cohort overlap.** The extent to which different GWAS rely on the same source
  databases was quantified and related to the strength of the diagnostic methods
  behind those databases, in order to estimate the risk of shared individuals
  between studies.
- **Genetic data heterogeneity.** An R pipeline standardises the genotyping
  platform manifests into a common format and crosses the resulting marker sets,
  measuring SNP overlap both between platforms and between studies.

## What this repository contains

- The full **R pipeline** for manifest standardisation and overlap analysis, from
  the raw platform files to the supplementary tables.
- The **Python notebook** used for the exploratory analysis and for the figures
  presented in the report.
- The **curated metadata table** and supplementary workbooks, provided for
  consultation.
- The **final figures** included in the report.

This repository accompanies a curricular internship report and supports a
methodological review of Alzheimer's disease GWAS currently in preparation by the
same team.

## Repository structure
alzheimer-gwas-curation/

├── R/                                     # Standardisation and analysis pipeline (R)

│   ├── pipeline.R                         # Master script, runs the full pipeline

│   ├── funcoes_padronizacao.R             # Manifest standardisation functions

│   ├── preparar_gwas.R                    # Cleans and structures the GWAS table

│   ├── analise_gwas.R                     # Union/intersection of SNPs per GWAS

│   ├── analise_gwas_cruzamento.R          # Pairwise SNP overlap between GWAS

│   ├── analise_cruzamento_platf.R         # Pairwise SNP overlap between platforms

│   ├── generate_supplementary_excel.R     # Builds the supplementary workbook

│   ├── heatmap_gwas.R                      # GWAS similarity heatmaps (run separately)

│   └── README.md                          # Notes on running the R pipeline

├── notebooks/                             # Exploratory analysis and figures (Python)

│   ├── graphs.ipynb                       # Cohort, diagnostic and SNP coverage figures

│   └── README.md                          # Notes on running the notebook

├── data/

│   └── raw/                               # Input data (see Data section)

├── supplementary/                         # Curated workbooks for consultation

├── figures/                               # Final figures used in the report

├── LICENSE.txt

└── README.md

## Data

The `data/raw/` folder contains the two input files used by the analysis,
`metadata_gwas.csv` and `metadata_snps.csv`. The hand-built metadata workbook in
`supplementary/` is the original source from which these files derive.

The genotyping platform manifests and annotation files (Illumina, Affymetrix and
Thermo Fisher) are **not included** in this repository, as they are subject to the
providers' own terms of use and cannot be redistributed. To reproduce the
standardisation step, they must be downloaded directly from the original sources:

- UCSC Genome Browser SNP array tables (https://genome.ucsc.edu/)
- Illumina product support documentation (https://support.illumina.com/)
- Thermo Fisher GeneChip array annotation files (https://www.thermofisher.com/)

The required filenames are listed in the header of `R/pipeline.R`.

## Requirements

The pipeline was developed in R. The following packages are required:

```r
install.packages(c("zoo", "DBI", "RSQLite", "openxlsx", "circlize", "RColorBrewer"))
```

`heatmap_gwas.R` also needs ComplexHeatmap, a Bioconductor package:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

The notebook was developed in Python 3 and uses `pandas`, `numpy` and
`matplotlib`.

## Running the pipeline

The full pipeline is orchestrated by `pipeline.R`, which runs every stage in
order, from standardising the platform manifests to generating the supplementary
workbook. It is computationally heavy, taking roughly 15 to 30 minutes, since it
processes all array manifests.

The pipeline runs from the `R/` folder. The scripts refer to one another and to
the intermediate files by plain filename, and read `metadata_gwas.csv` from
`data/raw/` through a relative path, so the working directory must be `R/`.

1. Download the platform manifests into the `R/` folder. The expected filenames
   are listed in the header of `pipeline.R`.
2. From R, set the working directory to `R/` and run:

```r
source("pipeline.R")
```

The input metadata is read from `data/raw/` automatically. The output is written
to the `R/` folder, namely the standardised marker files, the pairwise overlap
matrices and the supplementary Excel workbook.

## Figures

The figures in the report come from two sources. The Python notebook in
`notebooks/` produces the cohort, diagnostic and SNP coverage figures, writing its
output to `created_figuras/` and `created_tables/`. `R/heatmap_gwas.R` produces
the GWAS similarity heatmaps from the pairwise Jaccard matrices. The final figures
selected for the report are kept separately in `figures/`. See the README in each
folder for details.

## Author

Maria Leonor Soares Carvalho — curricular internship in Bioinformatics, i3S /
Faculty of Sciences, University of Porto (2026).

## License

This project is released under the MIT License. See `LICENSE.txt` for details.
