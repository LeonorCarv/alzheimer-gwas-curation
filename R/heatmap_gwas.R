# ============================================================================
# heatmap_gwas.R
# 
# Gera os heatmaps da semelhanca entre pares de GWAS, a partir das
# matrizes de Jaccard (uniao e intersecao) produzidas pelo analise_gwas_cruzamento.R.
# Usa o ComplexHeatmap, com agrupamento hierarquico e paleta Blues.
# ============================================================================

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

# ----------------------------------------------------------------------------
# 1. Ler as matrizes
# ----------------------------------------------------------------------------

m_union        <- as.matrix(read.table("matriz_gwas_uniao_jaccard.tsv",
                                       sep = "\t", header = TRUE,
                                       row.names = 1,
                                       check.names = FALSE))
m_intersection <- as.matrix(read.table("matriz_gwas_intersecao_jaccard.tsv",
                                       sep = "\t", header = TRUE,
                                       row.names = 1,
                                       check.names = FALSE))

cat("Union matrix:       ", nrow(m_union), "x", ncol(m_union), "\n")
cat("Intersection matrix:", nrow(m_intersection), "x", ncol(m_intersection), "\n")

# ----------------------------------------------------------------------------
# 2. Função que gera um heatmap a partir da matriz de Jaccard
# A distancia para o agrupamento é 1 - Jaccard, de modo a aproximar os
# estudos mais semelhantes. A matriz é apresentada em percentagem.
# ----------------------------------------------------------------------------

gerar_heatmap <- function(matriz, titulo, ficheiro_pdf, ficheiro_png) {
  
  # Converter em percentagem
  matriz_pct <- matriz * 100
  
  # Distância matriz e cluster
  matriz_dist <- as.dist(1 - matriz)
  cluster_obj <- hclust(matriz_dist, method = "complete")
  
  cores <- colorRampPalette(brewer.pal(9, "Blues"))(100)
  col_fun <- colorRamp2(seq(0, 100, length.out = 100), cores)
  
  # Heatmap
  ht <- Heatmap(
    matriz_pct,
    name             = "Similarity",
    col              = col_fun,
    cluster_rows     = cluster_obj,
    cluster_columns  = cluster_obj,
    show_row_names   = TRUE,
    show_column_names = TRUE,
    row_names_gp     = gpar(fontsize = 9),
    column_names_gp  = gpar(fontsize = 9),
    column_names_rot = 45,
    border           = FALSE,
    width            = unit(14 * ncol(matriz_pct), "pt"),
    height           = unit(14 * nrow(matriz_pct), "pt"),
    heatmap_legend_param = list(
      title           = "Similarity (based on Jaccard index, %)",
      title_position  = "leftcenter-rot",
      title_gp        = gpar(fontsize = 10, fontface = "plain"),
      labels_gp       = gpar(fontsize = 9),
      legend_height   = unit(6, "cm"),
      at              = c(0, 20, 40, 60, 80, 100),
      labels          = c("0", "20", "40", "60", "80", "100")
    )
  )
  
  # Guardar em formato em pdf com o título centrado
  pdf(ficheiro_pdf, width = 8, height = 7.3)
  draw(ht,
       heatmap_legend_side = "right",
       column_title        = titulo,
       column_title_gp     = gpar(fontsize = 13, fontface = "bold"),
       column_title_side   = "top")
  dev.off()
  
  # Guardar em formato png
  png(ficheiro_png, width = 8, height = 7.3, units = "in", res = 300)
  draw(ht,
       heatmap_legend_side = "right",
       column_title        = titulo,
       column_title_gp     = gpar(fontsize = 13, fontface = "bold"),
       column_title_side   = "top")
  dev.off()
  
  cat("Gerado:", ficheiro_pdf, "e", ficheiro_png, "\n")
}

# ----------------------------------------------------------------------------
# 3. Gerar os dois heatmaps
# ----------------------------------------------------------------------------

# Heatmap baseado na união
gerar_heatmap(
  matriz       = m_union,
  titulo       = "Union of SNP sets per GWAS\n",
  ficheiro_pdf = "heatmap_gwas_uniao.pdf",
  ficheiro_png = "heatmap_gwas_uniao.png"
)

# Heatmap baseado na interseção
gerar_heatmap(
  matriz       = m_intersection,
  titulo       = "Intersection of SNP sets per GWAS\n",
  ficheiro_pdf = "heatmap_gwas_intersecao.pdf",
  ficheiro_png = "heatmap_gwas_intersecao.png"
)

cat("\Ficheiros guardados:\n")
cat("  heatmap_gwas_uniao.pdf / .png\n")
cat("  heatmap_gwas_intersecao.pdf / .png\n")