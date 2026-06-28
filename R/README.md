# R scripts

Processing of the genotyping platforms and analysis of SNP overlap between
platforms and between GWAS, plus the GWAS similarity heatmaps. The scripts run in
a chain, each producing files that the next one reads.

## Execution order

`pipeline.R` runs the full analysis at once, in the correct order, and is the
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

## Figures

`heatmap_gwas.R` is run separately, after the pipeline. It reads the pairwise
Jaccard matrices produced in step 5 and draws the union and intersection
similarity heatmaps, with hierarchical clustering, saving them as PDF and PNG.

## Input data

The platform manifests (UCSC, Illumina, Axiom) and the GWAS metadata table
(`metadata_gwas.csv`, read from `data/raw/`). The manifests are not included in
the repository because they are large files from external sources. The header of
`pipeline.R` lists which ones they are and where to obtain them.

## How to run

The scripts refer to one another and to the intermediate files by plain filename,
and read `metadata_gwas.csv` from `data/raw/` through a relative path, so the
working directory must be `R/`. Place the downloaded manifests in `R/`, set the
working directory to `R/`, and run:

```r
source("pipeline.R")
```

To produce the heatmaps afterwards, with the same working directory, run:

```r
source("heatmap_gwas.R")
```

## Dependencies

Most of the pipeline runs on CRAN packages:

```r
install.packages(c("zoo", "DBI", "RSQLite", "openxlsx", "circlize", "RColorBrewer"))
```

`heatmap_gwas.R` also needs ComplexHeatmap, which is a Bioconductor package:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

The `grid` package ships with R and needs no installation.
