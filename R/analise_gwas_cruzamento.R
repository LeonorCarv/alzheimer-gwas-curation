# ============================================================================
# analise_gwas_cruzamento.R
#
# Calcula a matriz de SNPs partilhados entre GWAS (pairwise comparison).
#
# Gera 4 matrizes:
#   1. SNPs absolutos (uniao por GWAS)
#   2. Indice de Jaccard (uniao por GWAS)
#   3. SNPs absolutos (intersecao por GWAS)
#   4. Indice de Jaccard (intersecao por GWAS)
#
# Cada matriz e simetrica:
#   - Diagonal: total de SNPs do proprio GWAS
#   - Off-diagonal: SNPs partilhados ou indice de Jaccard
#
# Apenas inclui GWAS com pelo menos 1 subamostra valida.
# ============================================================================

# ----------------------------------------------------------------------------
# 1. Ler os conjuntos de SNPs por GWAS
# ----------------------------------------------------------------------------

cat("A ler conjuntos de SNPs por GWAS...\n")
conjuntos <- readRDS("conjuntos_snps_por_gwas.rds")

# Ler tambem o sumario por GWAS para ter autor e ano
sumario <- read.table("gwas_snps_sumario.tsv",
                      sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE,
                      na.strings = "NA",
                      quote = "", comment.char = "")

# Filtrar GWAS com pelo menos 1 subamostra valida
ids_validos <- sumario$ID[sumario$n_subamostras_com_plataforma > 0]
cat("GWAS com dados validos:", length(ids_validos), "\n")
cat("GWAS excluidos da matriz:", nrow(sumario) - length(ids_validos), "\n\n")

# Construir labels "Author Year" para os GWAS validos
sumario_validos <- sumario[sumario$ID %in% ids_validos, ]
sumario_validos <- sumario_validos[order(sumario_validos$ano,
                                         sumario_validos$autor), ]
labels_gwas <- paste(sumario_validos$autor, sumario_validos$ano)
ids_gwas_ordenados <- sumario_validos$ID

# Resolver duplicados de label (caso de Lambert 2013 duas vezes)
if (any(duplicated(labels_gwas))) {
  cat("AVISO: labels duplicados detectados. A acrescentar sufixo numerico.\n")
  labels_gwas <- make.unique(labels_gwas, sep = "_")
}

# ----------------------------------------------------------------------------
# 2. Construir listas de SNPs (uniao e intersecao) por GWAS
# ----------------------------------------------------------------------------

snps_uniao_por_gwas      <- list()
snps_intersecao_por_gwas <- list()

for (idx in seq_along(ids_gwas_ordenados)) {
  id <- ids_gwas_ordenados[idx]
  label <- labels_gwas[idx]
  snps_uniao_por_gwas[[label]]      <- conjuntos[[id]]$uniao
  snps_intersecao_por_gwas[[label]] <- conjuntos[[id]]$intersecao
}

# ----------------------------------------------------------------------------
# 3. Funcao para calcular matrizes (absoluta e Jaccard) a partir de uma lista
# ----------------------------------------------------------------------------

calcular_matrizes <- function(lista_snps, nome_conjunto) {
  n <- length(lista_snps)
  labels <- names(lista_snps)
  
  cat(sprintf("\nA calcular matrizes para '%s' (%d GWAS)...\n",
              nome_conjunto, n))
  
  matriz_absoluta <- matrix(0L, nrow = n, ncol = n,
                            dimnames = list(labels, labels))
  matriz_jaccard  <- matrix(0,  nrow = n, ncol = n,
                            dimnames = list(labels, labels))
  
  for (i in seq_len(n)) {
    snps_i <- lista_snps[[i]]
    n_i <- length(snps_i)
    
    for (j in seq_len(n)) {
      if (i == j) {
        matriz_absoluta[i, j] <- n_i
        matriz_jaccard[i, j]  <- 1
      } else if (i < j) {
        snps_j <- lista_snps[[j]]
        n_j <- length(snps_j)
        
        partilhados <- length(intersect(snps_i, snps_j))
        uniao_ij    <- n_i + n_j - partilhados
        
        if (uniao_ij == 0) {
          jaccard <- 0
        } else {
          jaccard <- partilhados / uniao_ij
        }
        
        matriz_absoluta[i, j] <- partilhados
        matriz_absoluta[j, i] <- partilhados
        matriz_jaccard[i, j]  <- jaccard
        matriz_jaccard[j, i]  <- jaccard
      }
    }
    cat(sprintf("  %d/%d - %s\n", i, n, labels[i]))
  }
  
  return(list(absoluta = matriz_absoluta,
              jaccard  = round(matriz_jaccard, 4)))
}

# ----------------------------------------------------------------------------
# 4. Calcular as 4 matrizes
# ----------------------------------------------------------------------------

matrizes_uniao      <- calcular_matrizes(snps_uniao_por_gwas, "UNIAO")
matrizes_intersecao <- calcular_matrizes(snps_intersecao_por_gwas, "INTERSECAO")

# ----------------------------------------------------------------------------
# 5. Gravar como TSV
# ----------------------------------------------------------------------------

cat("\nA gravar ficheiros TSV...\n")

gravar_matriz <- function(m, nome) {
  m_df <- cbind(GWAS = rownames(m), m)
  write.table(m_df, nome, sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat("  Guardado:", nome, "\n")
}

gravar_matriz(matrizes_uniao$absoluta,
              "matriz_gwas_uniao_absoluta.tsv")
gravar_matriz(matrizes_uniao$jaccard,
              "matriz_gwas_uniao_jaccard.tsv")
gravar_matriz(matrizes_intersecao$absoluta,
              "matriz_gwas_intersecao_absoluta.tsv")
gravar_matriz(matrizes_intersecao$jaccard,
              "matriz_gwas_intersecao_jaccard.tsv")

# ----------------------------------------------------------------------------
# 6. Estatistica resumida
# ----------------------------------------------------------------------------

cat("\n=== Estatistica resumida ===\n")

# Para a matriz Jaccard (uniao), top 10 pares mais similares
m <- matrizes_uniao$jaccard
diag(m) <- NA
indices <- which(upper.tri(m), arr.ind = TRUE)
pares <- data.frame(
  GWAS_A = rownames(m)[indices[, 1]],
  GWAS_B = colnames(m)[indices[, 2]],
  Jaccard = m[indices],
  SNPs_partilhados = matrizes_uniao$absoluta[indices],
  stringsAsFactors = FALSE
)
pares <- pares[order(-pares$Jaccard), ]

cat("\nTop 10 pares de GWAS mais similares (uniao, Jaccard):\n")
print(head(pares, 10))

cat("\nTop 10 pares de GWAS menos similares (uniao, Jaccard maior que zero):\n")
pares_nao_zero <- pares[pares$Jaccard > 0, ]
print(tail(pares_nao_zero, 10))

# Guardar o ranking de pares
write.table(pares, "ranking_pares_gwas_uniao.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\nGuardado: ranking_pares_gwas_uniao.tsv\n")

# Salvar objetos R para usar depois no Excel
saveRDS(list(matrizes_uniao      = matrizes_uniao,
             matrizes_intersecao = matrizes_intersecao,
             labels_gwas         = labels_gwas,
             ids_gwas            = ids_gwas_ordenados,
             ranking_pares       = pares),
        "matrizes_gwas.rds")
cat("Guardado: matrizes_gwas.rds (para uso no Excel)\n")

cat("\nAnalise de cruzamento entre GWAS concluida.\n")