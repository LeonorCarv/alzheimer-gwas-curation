# Notebooks

Visualizations for the Alzheimer's disease GWAS curation project. The notebook
starts from the study metadata tables and produces the figures and table used in
the analysis.

## Contents

- `graphs.ipynb` produces, in this order, the cohort-to-GWAS connection graphs
  (cases and controls), the diagnostic-method-to-GWAS connection graphs (cases
  and controls), the heatmap of shared individuals between studies, the table of
  cohorts by diagnostic method, and the stacked bar chart of SNP coverage per
  GWAS.

## Input data

The notebook reads two files from `data/raw/`, both semicolon-separated and
cp1252-encoded:

- `metadata_gwas.csv`, the GWAS metadata table
- `metadata_snps.csv`, the per-study SNP table

The paths in the reading cells are relative (`../data/raw/`), assuming the
notebook is run from the `notebooks/` folder.

## Generated files

When run, the notebook creates the output folders automatically:

- `created_figuras/`, the figures as 300 dpi PNGs, organised into subfolders by
  chart type (`sobrep_BD_GWAS`, `sobrep_Diag_GWAS`, `snpsVSgwas`)
- `created_tables/`, the table of studies by platform, as CSV and HTML

## How to run

Open the notebook in Jupyter and run the cells in order (Run menu, Run All).

## Dependencies

Python 3 with *pandas*, *NumPy* and *matplotlib*. Everything else is from the
standard library.
