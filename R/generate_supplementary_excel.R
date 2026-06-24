# ============================================================================
# generate_supplementary_excel.R
#
# Generates an Excel file with supplementary tables for the review of
# genotyping platforms used in GWAS.
#
# Structure of the file:
#   Sheet 1: README
#   Sheet 2: Included platforms summary
#   Sheet 3: Excluded platforms
#   Sheet 4: Methodological decisions
#   Sheet 5 - GWAS SNP analysis
#   Sheet 6 - GWAS overlap matrix (union, absolute)
#   Sheet 7 - GWAS overlap matrix (union, Jaccard)
#   Sheet 8 - GWAS overlap matrix (intersection, absolute)
#   Sheet 9 - GWAS overlap matrix (intersection, Jaccard)
#   Sheet 10 - GWAS pairs ranking

# Data for included platforms is automatically collected from padrao_*.tsv
# files present in the working directory.
# ============================================================================

library(openxlsx)

# ----------------------------------------------------------------------------
# 1. Collect data of included platforms from padrao files
# ----------------------------------------------------------------------------

# Metadata map per platform
metadata <- list(
  "padrao_Affy250K-Nsp.tsv" = list(
    platform   = "Affymetrix GeneChip Human Mapping 250K Nsp Array",
    manufacturer = "Affymetrix",
    source     = "UCSC (snpArrayAffy250Nsp)",
    version    = "snpArrayAffy250Nsp (UCSC)",
    notes      = ""
  ),
  "padrao_Affy250K-Sty.tsv" = list(
    platform   = "Affymetrix GeneChip Human Mapping 250K Sty Array",
    manufacturer = "Affymetrix",
    source     = "UCSC (snpArrayAffy250Sty)",
    version    = "snpArrayAffy250Sty (UCSC)",
    notes      = ""
  ),
  "padrao_AffySNP6.tsv" = list(
    platform   = "Applied Biosystems Human Genome-Wide SNP Array 6.0",
    manufacturer = "Affymetrix / Thermo Fisher",
    source     = "UCSC (snpArrayAffy6)",
    version    = "snpArrayAffy6 (UCSC)",
    notes      = "1,102 SNP_A markers originally annotated as 'unknown' included via positional matching"
  ),
  "padrao_Axiom-CEU.tsv" = list(
    platform   = "Applied Biosystems Axiom Genome-Wide CEU Array",
    manufacturer = "Thermo Fisher",
    source     = "Thermo Fisher (Axiom CEU SQLite)",
    version    = "Axiom_GW_Hu_SNP-na35-annot-db (dbSNP 142)",
    notes      = "Annotated against dbSNP 142 (Oct 2014), distinct from UCSC source (dbSNP 132)"
  ),
  "padrao_GSA-v1.tsv" = list(
    platform   = "Illumina Infinium Global Screening Array v1.0",
    manufacturer = "Illumina",
    source     = "Illumina (CSV manifest)",
    version    = "GSA-24v1-0_C1 (Jan 2018)",
    notes      = ""
  ),
  "padrao_Human610-Quad.tsv" = list(
    platform   = "Illumina Infinium Human610-Quad",
    manufacturer = "Illumina",
    source     = "Illumina (BED file)",
    version    = "Human610-Quadv1_H",
    notes      = "BED file lacks allele column; indels could not be filtered (documented limitation)"
  ),
  "padrao_Human660W-Quad.tsv" = list(
    platform   = "Illumina Infinium Human660W-Quad",
    manufacturer = "Illumina",
    source     = "UCSC (snpArrayIlluminaHuman660W_Quad)",
    version    = "snpArrayIlluminaHuman660W_Quad (UCSC)",
    notes      = ""
  ),
  "padrao_HumanCoreExome-12-v1-0.tsv" = list(
    platform   = "Illumina Infinium HumanCoreExome-12 v1.0",
    manufacturer = "Illumina",
    source     = "Illumina (CSV manifest + auxiliary file)",
    version    = "HumanCoreExome-12-v1-0-D (Apr 2014)",
    notes      = "Auxiliary file Name to RsID mapping applied to 254,432 markers"
  ),
  "padrao_HumanCoreExome-12-v1-1.tsv" = list(
    platform   = "Illumina Infinium HumanCoreExome-12 v1.1",
    manufacturer = "Illumina",
    source     = "Illumina (CSV manifest + auxiliary file)",
    version    = "HumanCoreExome-12-v1-1-C",
    notes      = "Auxiliary file Name to RsID mapping applied to 258,723 markers"
  ),
  "padrao_HumanExome-12-v1-0.tsv" = list(
    platform   = "Illumina Infinium HumanExome v1.0",
    manufacturer = "Illumina",
    source     = "Illumina (CSV manifest)",
    version    = "HumanExome-12-v1-0-B (Apr 2014)",
    notes      = "Also used to represent 'HumanExome' without version in ambiguous GWAS (ref. Grove et al. 2013)"
  ),
  "padrao_HumanExome-12-v1-1.tsv" = list(
    platform   = "Illumina Infinium HumanExome v1.1",
    manufacturer = "Illumina",
    source     = "Illumina (CSV manifest)",
    version    = "HumanExome-12-v1-1-B (Apr 2014)",
    notes      = ""
  ),
  "padrao_HumanHap300.tsv" = list(
    platform   = "Illumina Infinium HumanHap300",
    manufacturer = "Illumina",
    source     = "UCSC (snpArrayIllumina300)",
    version    = "snpArrayIllumina300 (UCSC)",
    notes      = ""
  ),
  "padrao_HumanHap550.tsv" = list(
    platform   = "Illumina Infinium HumanHap550",
    manufacturer = "Illumina",
    source     = "UCSC (snpArrayIllumina550)",
    version    = "snpArrayIllumina550 (UCSC)",
    notes      = ""
  ),
  "padrao_HumanHap650Y.tsv" = list(
    platform   = "Illumina Infinium HumanHap650Y",
    manufacturer = "Illumina",
    source     = "UCSC (snpArrayIllumina650)",
    version    = "snpArrayIllumina650 (UCSC)",
    notes      = ""
  ),
  "padrao_HumanOmni1-Quad.tsv" = list(
    platform   = "Illumina Infinium HumanOmni1-Quad",
    manufacturer = "Illumina",
    source     = "UCSC (snpArrayIlluminaHumanOmni1_Quad)",
    version    = "snpArrayIlluminaHumanOmni1_Quad (UCSC)",
    notes      = ""
  ),
  "padrao_HumanOmni2.5-Quad.tsv" = list(
    platform   = "Illumina Infinium HumanOmni2.5-Quad",
    manufacturer = "Illumina",
    source     = "Illumina (CSV manifest)",
    version    = "HumanOmni2.5-4v1_H (Apr 2011)",
    notes      = "Original build '37.1' normalised to '37' (hg19)"
  ),
  "padrao_HumanOmniExpress-12-v1-0.tsv" = list(
    platform   = "Illumina Infinium HumanOmniExpress-12 v1.0",
    manufacturer = "Illumina",
    source     = "Illumina (CSV manifest)",
    version    = "HumanOmniExpress-12-v1-0-K (Apr 2014)",
    notes      = ""
  ),
  "padrao_HumanOmniExpress-24-v1-0.tsv" = list(
    platform   = "Illumina Infinium HumanOmniExpress-24 v1.0",
    manufacturer = "Illumina",
    source     = "Illumina (CSV manifest)",
    version    = "HumanOmniExpress-24-v1-0-B (Apr 2014)",
    notes      = ""
  ),
  "padrao_HumanOmniExpress-24-v1-1.tsv" = list(
    platform   = "Illumina Infinium HumanOmniExpress-24 v1.1",
    manufacturer = "Illumina",
    source     = "Illumina (Annotation File)",
    version    = "HumanOmniExpress-24v1-1_A (Nov 2014)",
    notes      = "CSV manifest not available; Gene Annotation File used"
  ),
  "padrao_HumanOmniExpressExome-8-v1-2.tsv" = list(
    platform   = "Illumina Infinium HumanOmniExpressExome-8 v1.2",
    manufacturer = "Illumina",
    source     = "Illumina (CSV manifest)",
    version    = "HumanOmniExpressExome-8-v1-2-B (Apr 2014)",
    notes      = "Hybrid array (OmniExpress + HumanExome content)"
  ),
  "padrao_Illumina1M.tsv" = list(
    platform   = "Illumina Infinium Human1M-Duo",
    manufacturer = "Illumina",
    source     = "UCSC (snpArrayIllumina1M)",
    version    = "snpArrayIllumina1M (UCSC)",
    notes      = ""
  ),
  "padrao_MEGA.tsv" = list(
    platform   = "Illumina Infinium Multi-Ethnic Global-8 (MEGA)",
    manufacturer = "Illumina",
    source     = "Illumina (CSV manifest)",
    version    = "Multi-EthnicGlobal_D1 (Feb 2017)",
    notes      = "Also used to represent MEGAEX (MEGAEX manifest not publicly available)"
  ),
  "padrao_PsychArray.tsv" = list(
    platform   = "Illumina Infinium PsychArray",
    manufacturer = "Illumina",
    source     = "Illumina (Annotation File + auxiliary file)",
    version    = "PsychArray_A (Mar 2014)",
    notes      = "CSV manifest not available; Gene Annotation File and auxiliary file used (279,111 mappings applied)"
  )
)

# Collect counts from each padrao file
cat("Reading padrao files...\n")
included_summary <- data.frame(
  N                   = integer(),
  Platform            = character(),
  Manufacturer        = character(),
  Primary_source      = character(),
  Manifest_version    = character(),
  Build               = character(),
  Total_markers       = integer(),
  Type_rs             = integer(),
  Type_positional     = integer(),
  Percentage_rs       = numeric(),
  Notes               = character(),
  stringsAsFactors    = FALSE
)

i <- 1
for (file in names(metadata)) {
  if (!file.exists(file)) {
    cat("WARNING: file not found:", file, "\n")
    next
  }
  
  data <- read.table(file, sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE,
                     na.strings = "NA", quote = "", comment.char = "")
  
  n_total <- nrow(data)
  n_rs    <- sum(data$tipo_id == "rs")
  n_pos   <- sum(data$tipo_id == "posicao")
  meta    <- metadata[[file]]
  
  included_summary <- rbind(included_summary, data.frame(
    N                = i,
    Platform         = meta$platform,
    Manufacturer     = meta$manufacturer,
    Primary_source   = meta$source,
    Manifest_version = meta$version,
    Build            = "hg19 (GRCh37)",
    Total_markers    = n_total,
    Type_rs          = n_rs,
    Type_positional  = n_pos,
    Percentage_rs    = round(n_rs / n_total * 100, 2),
    Notes            = meta$notes,
    stringsAsFactors = FALSE
  ))
  i <- i + 1
}

# Add MEGAEX row (represented by MEGA)
mega_row <- included_summary[included_summary$Platform == 
                               "Illumina Infinium Multi-Ethnic Global-8 (MEGA)", ]
if (nrow(mega_row) > 0) {
  megaex_row <- mega_row
  megaex_row$N <- max(included_summary$N) + 1
  megaex_row$Platform <- "Illumina Infinium Expanded Multi-Ethnic Genotyping Array (MEGAEX)"
  megaex_row$Primary_source <- "Approximated by MEGA"
  megaex_row$Manifest_version <- "Multi-EthnicGlobal_D1 (Feb 2017) used by approximation"
  megaex_row$Notes <- "MEGAEX manifest not publicly available; represented by MEGA manifest"
  included_summary <- rbind(included_summary, megaex_row)
}

cat("Included platforms:", nrow(included_summary), "\n")

# ----------------------------------------------------------------------------
# 2. Excluded platforms table
# ----------------------------------------------------------------------------

excluded <- data.frame(
  N        = 1:16,
  Platform = c(
    "Illumina Infinium HumanCNV370-Duo",
    "Illumina Infinium Neurochip",
    "Applied Biosystems Axiom Spain Biobank Array",
    "Applied Biosystems UK BiLEVE Axiom Array",
    "Applied Biosystems UK Biobank Axiom Array",
    "Illumina GoldenGate Custom Genotyping",
    "Illumina iSelect Custom Genotyping BeadChip",
    "Illumina Infinium Metabochip",
    "Illumina Infinium UM HUNT Biobank v1.0",
    "deCODE genetics Custom Axiom Array v1",
    "Agena Bioscience MassARRAY System",
    "LGC KASP Genotyping Assay",
    "Sequenom MassARRAY System",
    "TaqMan Genotyping Assay",
    "Illumina HiSeqX Ten",
    "Illumina NovaSeq 6000"
  ),
  Category = c(
    "Catalog array",
    "Catalog array / Consortium",
    "Customised",
    "Customised",
    "Customised",
    "Customised",
    "Customised",
    "Semi-custom",
    "Customised / Consortium",
    "Customised",
    "Non-array chemistry",
    "Non-array chemistry",
    "Non-array chemistry",
    "Non-array chemistry",
    "Sequencing",
    "Sequencing"
  ),
  Exclusion_reason = c(
    "Manifest available only in proprietary binary BPM format; no open-source reader accessible",
    "Consortium-designed array (IPDGC); no publicly accessible manifest",
    "Custom; no publicly accessible manifest",
    "Custom; no publicly accessible manifest",
    "Custom; no publicly accessible manifest",
    "Generic technology for custom arrays; no single manifest exists",
    "Generic technology for custom arrays; no single manifest exists",
    "Semi-custom (Illumina + Metabochip Consortium); manifest not public",
    "Consortium-designed array; manifest not public",
    "Custom (deCODE Genetics); manifest private",
    "Non-array chemistry; hardware-level platform comparison not applicable",
    "Non-array chemistry; hardware-level platform comparison not applicable",
    "Non-array chemistry; hardware-level platform comparison not applicable",
    "Non-array chemistry; hardware-level platform comparison not applicable",
    "Sequencing platform; not applicable to array content comparison",
    "Sequencing platform; not applicable to array content comparison"
  ),
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------------------------
# 3. Methodological decisions table
# ----------------------------------------------------------------------------

decisions <- data.frame(
  ID = 1:18,
  Category = c(
    "Inclusion",
    "Inclusion",
    "Inclusion",
    "Inclusion",
    "Inclusion",
    "Data source",
    "Data source",
    "Reference genome",
    "Reference genome",
    "Content filtering",
    "Content filtering",
    "Content filtering",
    "Content filtering",
    "Content filtering",
    "Processing",
    "Processing",
    "Processing",
    "Cross-platform comparison"
  ),
  Decision = c(
    "Inclusion of catalog arrays only (commercial off-the-shelf arrays)",
    "Exclusion of customised platforms (Axiom Spain, UK BiLEVE, UK Biobank, GoldenGate, iSelect, Metabochip, Neurochip, HUNT, deCODE)",
    "Exclusion of non-array chemistry (MassARRAY, KASP, TaqMan) and sequencing platforms (HiSeqX, NovaSeq)",
    "Version ambiguity resolution: OmniExpress-24 without specified version assumed as v1.3 for GWAS published 2019-2022",
    "Documented substitution: MEGAEX represented by MEGA manifest, given the unavailability of a public MEGAEX manifest",
    "Source hierarchy: UCSC first (snpArray tables), Illumina/Thermo Fisher for arrays not available in UCSC",
    "Recognition of dbSNP version differences across sources (UCSC: build 131-132; Axiom: build 142; Illumina: variable). dbSNP version harmonisation deferred to a subsequent step",
    "Adoption of hg19 (GRCh37) as the single reference assembly of the project",
    "One-based coordinate convention adopted across the project. Conversion from zero-based (UCSC, BED) to one-based applied at processing time",
    "Restriction to SNPs in the strict sense: only variants with canonical alleles A, C, G, T accepted. Indels (D/I, I/D, hyphen) excluded",
    "Restriction to nuclear chromosomes (1-22, X, Y). Excluded: mitochondrial (chrM, 26 in Axiom), unplaced scaffolds (chrN_random, chrUn_), markers without mapping (INT_MAX in Axiom files)",
    "Exclusion of copy number variation (CNV) content. Identified by the 'cnvi' prefix in Illumina tables and by a separate table in Affymetrix arrays",
    "Inclusion of markers with proprietary identifiers (GA, VG, SNP_other, hCV, pgxUn, numeric, FY, DI, HPA, JK) via positional matching, when meeting valid-nuclear-SNP criteria",
    "AffySNP6 markers with rsId 'unknown' (1,102 markers) included via positional matching. The string 'unknown' was replaced by NA in the original-name column",
    "Construction of pseudo-observed-allele column for Axiom arrays (from Allele_A/Allele_B) and Illumina arrays (from SNP column in [A/G] format)",
    "Extraction of rs identifiers embedded in prefixes of the type 'exm-rs<number>' (case of HumanExome v1.0 and v1.1 and other arrays)",
    "Application of auxiliary mapping files Name to RsID only for non-rs Names with a single valid rs (1% threshold criterion). Excluded: harmonisation of old rs to new rs, discontinued RsID ('.'), multiple alternative rs identifiers",
    "Explicit limitation: comparison is performed at the level of the nominal content of each array (SNPs printed on probes), not at the level of SNPs effectively tested in each study (which would depend on imputation and quality control)"
  ),
  Justification = c(
    "Lack of public manifests for customised arrays; ad-hoc content does not allow methodologically robust comparison",
    "Same as above",
    "Platforms based on different chemistry do not allow hardware-level array comparison",
    "Dominant version during the publication period of the studies. Assumption documented as a limitation",
    "MEGA and MEGAEX belong to the same family with substantial content overlap; approximation with documented limitation",
    "Maximise internal consistency in annotation source, given that UCSC arrays share the same annotation date",
    "Inherent differences across data sources; harmonisation requires a dedicated step",
    "Most commonly used build during the publication period of the analysed GWAS and across the data sources used",
    "Convention used by Illumina, Thermo Fisher, dbSNP, and VCF format",
    "Ensure comparability across platforms; positional matching does not distinguish SNP from indel",
    "Focus on nuclear variants; scaffold and mitochondrial variants are residual for most arrays",
    "Such markers belong to the functional CNV programme (copy number variation), not the SNP programme",
    "Consistency with valid-nuclear-SNP criterion; preservation of traceability in the original-name column",
    "Real markers without dbSNP rs assigned at annotation time. Consistency with the general rule",
    "Allows uniform application of the canonical-SNP filter to all arrays",
    "Recovery of canonical rs identifiers hidden behind manufacturer proprietary prefixes",
    "Maximisation of rs coverage without introducing dbSNP version harmonisation (handled in a subsequent step)",
    "Imputation and QC are study-specific; the comparison of array content is the focus of this study"
  ),
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------------------------
# 4. README sheet with introduction
# ----------------------------------------------------------------------------

readme <- data.frame(
  Section = c(
    "Title",
    "Aim",
    "Period of analysed GWAS",
    "Reference genome",
    "Coordinate convention",
    "Contents of this file",
    "",
    "Sheet 2 - Included platforms summary",
    "Sheet 3 - Excluded platforms",
    "Sheet 4 - Methodological decisions",
    "Sheet 5 - GWAS SNP analysis",
    "Sheet 6 - GWAS_union_abs",
    "Sheet 7 - GWAS_union_jaccard",
    "Sheet 8 - GWAS_intersection_abs",
    "Sheet 9 - GWAS_intersection_jaccard",
    "Sheet 10 - GWAS pairs ranking",
    "",
    "Definitions",
    "  rs",
    "  positional",
    "  hg19 / GRCh37",
    "  Build 37.1",
    "  CNV (cnvi)",
    "  Indels",
    "  Canonical SNP",
    "",
    "Extraction date",
    "Author"
  ),
  Description = c(
    "Supplementary tables: comparison of genotyping platforms used in GWAS",
    "Identify SNPs shared and unique across genotyping platforms used in genome-wide association studies",
    "2019, 2021, and 2022",
    "hg19 (equivalent to GRCh37)",
    "one-based (the first base of the chromosome is numbered 1)",
    "This file contains detailed information about the genotyping platforms included in and excluded from the study, together with the methodological decisions made during processing.",
    "",
    "List of commercial platforms included in the study, with marker counts, manifest source, exact version, and relevant notes",
    "List of excluded platforms, with category and reason for exclusion",
    "List of methodological decisions made during processing, with justification for each",
    "Per-GWAS analysis of SNP coverage. For each GWAS, the table shows the number of subsamples, the number of platforms used (and which), and two SNP counts: the UNION (SNPs in at least one subsample) and the INTERSECTION (SNPs available in ALL subsamples with platform known, the biologically analysable set before imputation). Subsamples with platform unknown are excluded from both calculations. The Has_complete_data column flags GWAS where every subsample has a known platform; these are the most robust for cross-study comparison.",
    "Pairwise matrix of SNP overlap between GWAS, based on the UNION sets. Each cell shows the number of SNPs shared between two GWAS. Diagonal shows total SNPs per GWAS. Only GWAS with at least one valid subsample are included.",
    "Pairwise matrix of Jaccard similarity index between GWAS, based on the UNION sets. Jaccard index ranges from 0 (no overlap) to 1 (identical content). Calculated as |A intersection B| / |A union B|. Diagonal values are 1.",
    "Pairwise matrix of SNP overlap between GWAS, based on the INTERSECTION sets (SNPs present in ALL subsamples). Represents the realistic shared content for meta-analysis without imputation.",
    "Pairwise matrix of Jaccard similarity index between GWAS, based on the INTERSECTION sets. Same scale as Sheet 7 but using only SNPs analysable in all subsamples per GWAS.",
    "Ranking of all GWAS pairs ordered by similarity. Columns include both union and intersection metrics (shared SNPs and Jaccard) to allow flexible comparison. Use Excel filters to find the most similar or least similar pairs of studies.",
    "",
    "",
    "dbSNP variant identifier (of the form rs123456789)",
    "Marker without dbSNP rs assigned, identified by chromosomal position only",
    "Genome Reference Consortium Human Build 37, released in February 2009",
    "Build 37 variant with minor corrections. Considered equivalent to hg19 for most applications",
    "Copy number variation. Specific type of marker, excluded from this study",
    "Insertions and deletions. Marker type excluded from this study",
    "SNP in which both alleles are canonical nucleotides (A, C, G, or T). Distinguished from indels and multi-base variants",
    "",
    as.character(Sys.Date()),
    "[Author name to be filled in]"
  ),
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------------------------
# 5. Build the Excel file
# ----------------------------------------------------------------------------

wb <- createWorkbook()

# Styles
header_style <- createStyle(
  fontColour = "white",
  fgFill     = "#4F81BD",
  textDecoration = "bold",
  halign     = "center",
  valign     = "center",
  border     = "TopBottomLeftRight",
  wrapText   = TRUE
)

cell_style <- createStyle(
  valign   = "top",
  wrapText = TRUE,
  border   = "TopBottomLeftRight"
)

# Sheet 1: README
addWorksheet(wb, "README")
writeData(wb, "README", readme, startCol = 1, startRow = 1)
addStyle(wb, "README", header_style, rows = 1, cols = 1:2)
addStyle(wb, "README", cell_style, 
         rows = 2:(nrow(readme)+1), cols = 1:2,
         gridExpand = TRUE)
setColWidths(wb, "README", cols = 1, widths = 40)
setColWidths(wb, "README", cols = 2, widths = 100)
freezePane(wb, "README", firstRow = TRUE)

# Sheet 2: Included platforms
addWorksheet(wb, "Included_platforms")
writeData(wb, "Included_platforms", included_summary, startCol = 1, startRow = 1)
addStyle(wb, "Included_platforms", header_style,
         rows = 1, cols = 1:ncol(included_summary))
addStyle(wb, "Included_platforms", cell_style,
         rows = 2:(nrow(included_summary)+1), cols = 1:ncol(included_summary),
         gridExpand = TRUE)
setColWidths(wb, "Included_platforms", cols = 1, widths = 5)
setColWidths(wb, "Included_platforms", cols = 2, widths = 60)
setColWidths(wb, "Included_platforms", cols = 3, widths = 20)
setColWidths(wb, "Included_platforms", cols = 4, widths = 40)
setColWidths(wb, "Included_platforms", cols = 5, widths = 45)
setColWidths(wb, "Included_platforms", cols = 6, widths = 15)
setColWidths(wb, "Included_platforms", cols = 7:10, widths = 15)
setColWidths(wb, "Included_platforms", cols = 11, widths = 70)
freezePane(wb, "Included_platforms", firstRow = TRUE)
addFilter(wb, "Included_platforms", row = 1, cols = 1:ncol(included_summary))

# Sheet 3: Excluded platforms
addWorksheet(wb, "Excluded_platforms")
writeData(wb, "Excluded_platforms", excluded, startCol = 1, startRow = 1)
addStyle(wb, "Excluded_platforms", header_style,
         rows = 1, cols = 1:ncol(excluded))
addStyle(wb, "Excluded_platforms", cell_style,
         rows = 2:(nrow(excluded)+1), cols = 1:ncol(excluded),
         gridExpand = TRUE)
setColWidths(wb, "Excluded_platforms", cols = 1, widths = 5)
setColWidths(wb, "Excluded_platforms", cols = 2, widths = 50)
setColWidths(wb, "Excluded_platforms", cols = 3, widths = 25)
setColWidths(wb, "Excluded_platforms", cols = 4, widths = 90)
freezePane(wb, "Excluded_platforms", firstRow = TRUE)
addFilter(wb, "Excluded_platforms", row = 1, cols = 1:ncol(excluded))

# Sheet 4: Methodological decisions
addWorksheet(wb, "Methodological_decisions")
writeData(wb, "Methodological_decisions", decisions, startCol = 1, startRow = 1)
addStyle(wb, "Methodological_decisions", header_style,
         rows = 1, cols = 1:ncol(decisions))
addStyle(wb, "Methodological_decisions", cell_style,
         rows = 2:(nrow(decisions)+1), cols = 1:ncol(decisions),
         gridExpand = TRUE)
setColWidths(wb, "Methodological_decisions", cols = 1, widths = 5)
setColWidths(wb, "Methodological_decisions", cols = 2, widths = 25)
setColWidths(wb, "Methodological_decisions", cols = 3, widths = 80)
setColWidths(wb, "Methodological_decisions", cols = 4, widths = 80)
freezePane(wb, "Methodological_decisions", firstRow = TRUE)
addFilter(wb, "Methodological_decisions", row = 1, cols = 1:ncol(decisions))


# ----------------------------------------------------------------------------
# Sheet 5: GWAS - SNP analysis per study
# ----------------------------------------------------------------------------

cat("\nPreparing GWAS sheet...\n")

gwas_sumario <- read.table("gwas_snps_sumario.tsv",
                           sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE,
                           na.strings = "NA",
                           quote = "", comment.char = "")

# Calculate additional columns for readability
gwas_sumario$N_platforms_used <- sapply(
  strsplit(gwas_sumario$plataformas_unicas, ";\\s*"),
  function(x) length(x[x != ""])
)

gwas_sumario$Has_complete_data <- gwas_sumario$n_subamostras_NA == 0 &
  gwas_sumario$n_subamostras_com_plataforma > 0

# Build final table with Author and Year in separate columns
gwas_excel <- data.frame(
  First_author                  = gwas_sumario$autor,
  Year                          = gwas_sumario$ano,
  N_subsamples_total            = gwas_sumario$n_subamostras_total,
  N_subsamples_with_platform    = gwas_sumario$n_subamostras_com_plataforma,
  N_subsamples_unknown          = gwas_sumario$n_subamostras_NA,
  N_platforms_used              = gwas_sumario$N_platforms_used,
  Platforms_used                = gsub(";\\s*", ", ",
                                       gwas_sumario$plataformas_unicas),
  N_SNPs_union                  = gwas_sumario$n_snps_uniao,
  N_SNPs_intersection           = gwas_sumario$n_snps_intersecao,
  Has_complete_data             = gwas_sumario$Has_complete_data,
  stringsAsFactors              = FALSE
)

# Sort by year (oldest first) then alphabetically by author
gwas_excel <- gwas_excel[order(gwas_excel$Year, gwas_excel$First_author), ]

addWorksheet(wb, "GWAS_SNP_analysis")
writeData(wb, "GWAS_SNP_analysis", gwas_excel, startCol = 1, startRow = 1)
addStyle(wb, "GWAS_SNP_analysis", header_style,
         rows = 1, cols = 1:ncol(gwas_excel))
addStyle(wb, "GWAS_SNP_analysis", cell_style,
         rows = 2:(nrow(gwas_excel)+1), cols = 1:ncol(gwas_excel),
         gridExpand = TRUE)
setColWidths(wb, "GWAS_SNP_analysis", cols = 1, widths = 20)
setColWidths(wb, "GWAS_SNP_analysis", cols = 2, widths = 8)
setColWidths(wb, "GWAS_SNP_analysis", cols = 3:6, widths = 18)
setColWidths(wb, "GWAS_SNP_analysis", cols = 7, widths = 80)
setColWidths(wb, "GWAS_SNP_analysis", cols = 8:9, widths = 22)
setColWidths(wb, "GWAS_SNP_analysis", cols = 10, widths = 18)
freezePane(wb, "GWAS_SNP_analysis", firstRow = TRUE)
addFilter(wb, "GWAS_SNP_analysis", row = 1, cols = 1:ncol(gwas_excel))


# ----------------------------------------------------------------------------
# Sheet 6 to 9: GWAS x GWAS overlap matrices and ranking
# ----------------------------------------------------------------------------

cat("\nPreparing GWAS x GWAS matrices...\n")

matrizes_dados <- readRDS("matrizes_gwas.rds")

# Helper function to add a matrix sheet to the workbook
add_matrix_sheet <- function(wb, sheet_name, matrix_data, header_style, cell_style) {
  matrix_df <- cbind(GWAS = rownames(matrix_data), as.data.frame(matrix_data))
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, matrix_df, startCol = 1, startRow = 1)
  addStyle(wb, sheet_name, header_style,
           rows = 1, cols = 1:ncol(matrix_df))
  addStyle(wb, sheet_name, cell_style,
           rows = 2:(nrow(matrix_df)+1), cols = 1:ncol(matrix_df),
           gridExpand = TRUE)
  setColWidths(wb, sheet_name, cols = 1, widths = 25)
  setColWidths(wb, sheet_name, cols = 2:ncol(matrix_df), widths = 15)
  freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
}

# Sheet 6: Overlap matrix (union, absolute SNP counts)
add_matrix_sheet(wb, "GWAS_union_abs",
                 matrizes_dados$matrizes_uniao$absoluta,
                 header_style, cell_style)

# Sheet 7: Overlap matrix (union, Jaccard index)
add_matrix_sheet(wb, "GWAS_union_jaccard",
                 matrizes_dados$matrizes_uniao$jaccard,
                 header_style, cell_style)

# Sheet 8: Overlap matrix (intersection, absolute SNP counts)
add_matrix_sheet(wb, "GWAS_intersection_abs",
                 matrizes_dados$matrizes_intersecao$absoluta,
                 header_style, cell_style)

# Sheet 9: Overlap matrix (intersection, Jaccard index)
add_matrix_sheet(wb, "GWAS_intersection_jaccard",
                 matrizes_dados$matrizes_intersecao$jaccard,
                 header_style, cell_style)

# Sheet 10: Ranking of GWAS pairs (most and least similar)
cat("\nPreparing ranking of GWAS pairs...\n")

ranking <- matrizes_dados$ranking_pares
colnames(ranking) <- c("GWAS_A", "GWAS_B", "Jaccard_union",
                       "Shared_SNPs_union")

# Add Jaccard intersection and absolute intersection counts
m_int_jac <- matrizes_dados$matrizes_intersecao$jaccard
m_int_abs <- matrizes_dados$matrizes_intersecao$absoluta

ranking$Jaccard_intersection <- mapply(function(a, b) {
  if (a %in% rownames(m_int_jac) && b %in% colnames(m_int_jac)) {
    m_int_jac[a, b]
  } else NA
}, ranking$GWAS_A, ranking$GWAS_B)

ranking$Shared_SNPs_intersection <- mapply(function(a, b) {
  if (a %in% rownames(m_int_abs) && b %in% colnames(m_int_abs)) {
    m_int_abs[a, b]
  } else NA
}, ranking$GWAS_A, ranking$GWAS_B)

# Reorder columns for clarity
ranking <- ranking[, c("GWAS_A", "GWAS_B",
                       "Shared_SNPs_union", "Jaccard_union",
                       "Shared_SNPs_intersection", "Jaccard_intersection")]

addWorksheet(wb, "GWAS_pairs_ranking")
writeData(wb, "GWAS_pairs_ranking", ranking, startCol = 1, startRow = 1)
addStyle(wb, "GWAS_pairs_ranking", header_style,
         rows = 1, cols = 1:ncol(ranking))
addStyle(wb, "GWAS_pairs_ranking", cell_style,
         rows = 2:(nrow(ranking)+1), cols = 1:ncol(ranking),
         gridExpand = TRUE)
setColWidths(wb, "GWAS_pairs_ranking", cols = 1:2, widths = 25)
setColWidths(wb, "GWAS_pairs_ranking", cols = 3:6, widths = 22)
freezePane(wb, "GWAS_pairs_ranking", firstRow = TRUE)
addFilter(wb, "GWAS_pairs_ranking", row = 1, cols = 1:ncol(ranking))


# ----------------------------------------------------------------------------
# Save the file
# ----------------------------------------------------------------------------
output_file <- "Supplementary_tables_genotyping_platforms.xlsx"
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("\nExcel file created:", output_file, "\n")
cat("Sheets:\n")
cat("  1. README\n")
cat("  2. Included_platforms (", nrow(included_summary), "rows)\n")
cat("  3. Excluded_platforms (", nrow(excluded), "rows)\n")
cat("  4. Methodological_decisions (", nrow(decisions), "rows)\n")
cat("  5. GWAS_SNP_analysis (", nrow(gwas_excel), "rows)\n")
cat("  6. GWAS_union_abs (matrix",
    nrow(matrizes_dados$matrizes_uniao$absoluta), "x",
    ncol(matrizes_dados$matrizes_uniao$absoluta), ")\n")
cat("  7. GWAS_union_jaccard (same dimensions)\n")
cat("  8. GWAS_intersection_abs (same dimensions)\n")
cat("  9. GWAS_intersection_jaccard (same dimensions)\n")
cat("  10. GWAS_pairs_ranking (", nrow(ranking), "rows)\n")


