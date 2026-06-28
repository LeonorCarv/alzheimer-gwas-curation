# R scripts

Processing of the genotyping platforms and analysis of SNP overlap between
platforms and between GWAS. The scripts run in a chain, each producing files that
the next one reads.

## Execution order

`pipeline.R` runs everything at once, in the correct order, and is the
recommended entry point. Under the hood, the sequence is:

1. `funcoes_padronizacao.R` holds the functions that normalise the manifests from
   the various sources (UCSC, Illumina, Thermo Fisher Axiom) and merge them into a
   single master file.
2. `analise_cruzamento_platf.R` computes the SNP overlap between platforms.
3. `preparar_gwas.R` structures the GWAS table from the metadata CSV.
4. `analise_gwas.R` computes, for each GWAS, the SNPs by union and by intersection
   of the platforms used.
5. `analise_gwas_cruzamento.R` builds the pairwise SNP overlap matrices between
   GWAS.
6. `generate_supplementary_excel.R` gathers everything into a single Excel file
   with the supplementary tables.

## Input data

The platform manifests (UCSC, Illumina, Axiom) and the GWAS metadata table. The
manifests are not included in the repository because they are large files from
external sources. The header of `pipeline.R` lists which ones they are and where
to obtain them.

## How to run

From the project root, in R:

    source("R/pipeline.R")

The file paths are relative, so it is best to always run from the same folder.

## Dependencies

R with the *zoo*, *DBI*, *RSQLite* and *openxlsx* packages.
