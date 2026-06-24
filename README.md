# alzheimer-gwas-curation

Curation, standardisation and overlap analysis of metadata and genotyping platforms across forty Alzheimer's disease GWAS.

This repository contains the curation pipeline and analysis code developed during a curricular internship at i3S (Instituto de Investigação e Inovação
em Saúde), in support of a methodological review of genome-wide association studies (GWAS) in Alzheimer's disease. The work quantifies two forms of
overlap between studies, at the level of source cohorts and at the level of
genetic markers (SNPs), and characterises the heterogeneity of the genotyping platforms used.

## Repository structure
alzheimer-gwas-curation/

├── R/                  # Standardisation and analysis pipeline (R)
│   ├── pipeline.R                    # Master script, runs the full pipeline
│   ├── funcoes_padronizacao.R        # Manifest standardisation functions
│   ├── preparar_gwas.R               # Cleans and structures the GWAS table
│   ├── analise_gwas.R                # Union/intersection of SNPs per GWAS
│   ├── analise_gwas_cruzamento.R     # Pairwise SNP overlap between GWAS
│   ├── analise_cruzamento_platf.R    # Pairwise SNP overlap between platforms
│   └── generate_supplementary_excel.R # Builds the supplementary workbook
├── notebooks/          # Exploratory analysis and figures (Python)
├── data/
│   └── raw/            # Input data (see Data section)
├── supplementary/      # Curated workbooks for consultation
├── figures/            # Final figures used in the report
├── LICENSE.txt
└── README.md

## Data

The `data/raw/` folder contains the two input files used by the analysis,
`metadata_gwas.csv` and `metadata_snps.csv`. The hand-built metadata
workbook in `supplementary/` is the original source from which these files
derive.

The genotyping platform manifests and annotation files (Illumina,
Affymetrix and Thermo Fisher) are **not included** in this repository, as
they are subject to the providers' own terms of use and cannot be
redistributed. To reproduce the standardisation step, they must be
downloaded directly from the original sources:

- UCSC Genome Browser SNP array tables (https://genome.ucsc.edu/)
- Illumina product support documentation (https://support.illumina.com/)
- Thermo Fisher GeneChip array annotation files (https://www.thermofisher.com/)

The required filenames are listed in the header of `R/pipeline.R`.

## Requirements

The pipeline was developed in R. The following packages are required:

```r
install.packages(c("zoo", "DBI", "RSQLite", "openxlsx"))
```

The notebooks were developed in Python 3 and use `pandas`, `numpy` and
`matplotlib`.

## Running the pipeline

The full pipeline is orchestrated by `R/pipeline.R`, which runs every stage
in order, from standardising the platform manifests to generating the
supplementary workbook. Note that it is computationally heavy, taking
roughly 15 to 30 minutes, since it processes all array manifests.

1. Place the downloaded manifests and the input data in the working
   directory (the filenames expected are listed in the header of
   `pipeline.R`).
2. From R, with the working directory set to the project folder, run:

```r
source("R/pipeline.R")
```

The pipeline produces the standardised marker files, the pairwise overlap
matrices and the supplementary Excel workbook.

## Figures

The `notebooks/` folder contains the Python code used for the exploratory
analysis and for the figures presented in the report. The final versions of
those figures are available in `figures/`.

## Author

Maria Leonor Soares Carvalho — curricular internship in Bioinformatics,
i3S / Faculty of Sciences, University of Porto (2026).

## License

This project is released under the MIT License. See `LICENSE.txt` for
details.
