# ============================================================================
# analise_gwas_cruzamento.R
#
# Calcula a matriz de SNPs partilhados entre GWAS (pairwise comparison).
#
# Gera 4 matrizes:
#   1. SNPs partilhados em valor absolutos (união por GWAS)
#   2. Índice de Jaccard (uniao por GWAS)
#   3. SNPs partilhados em valor absolutos (interseção por GWAS)
#   4. Índice de Jaccard (intersecao por GWAS)
#
# Em cada matriz a diagonal guarda o total de SNPs do proprio GWAS e as
# restantes celulas o valor partilhado. So entram GWAS com pelo menos uma
# subamostra valida.
#
# Lê conjuntos_snps_por_gwas.rds e gwas_snps_sumario.tsv
# ============================================================================

# ----------------------------------------------------------------------------
# 1. Ler os conjuntos de SNPs por GWAS
# ----------------------------------------------------------------------------

conjuntos <- readRDS("conjuntos_snps_por_gwas.rds")

# Ler também o sumário por GWAS para ter autor e ano
sumario <- read.table("gwas_snps_sumario.tsv",
                      sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE,
                      na.strings = "NA",
                      quote = "", comment.char = "")

# Filtrar GWAS com pelo menos 1 subamostra válida
ids_validos <- sumario$ID[sumario$n_subamostras_com_plataforma > 0]
cat("GWAS com dados validos:", length(ids_validos), "\n")
cat("GWAS excluidos da matriz:", nrow(sumario) - length(ids_validos), "\n\n")

# Ordenar por ano e autor e montar os labels "Autor Ano"
sumario_validos <- sumario[sumario$ID %in% ids_validos, ]
sumario_validos <- sumario_validos[order(sumario_validos$ano,
                                         sumario_validos$autor), ]
labels_gwas <- paste(sumario_validos$autor, sumario_validos$ano)
ids_gwas_ordenados <- sumario_validos$ID

# Dois GWAS de Lambert 2013 partilham o mesmo label. Quando isso acontece,
# acrescenta-se um sufixo numerico para os distinguir nas matrizes.
if (any(duplicated(labels_gwas))) {
  cat("Labels duplicados detetados. Acrescentado sufixo numérico.\n")
  labels_gwas <- make.unique(labels_gwas, sep = "_")
}

# ----------------------------------------------------------------------------
# 2. Construir listas de SNPs (uniao e intersecao) por GWAS (indexados pelo label)
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
# 3. Função que devolve a matriz de partilhados em valor absoluto e a matriz 
# de Jaccard. A diagonal fica com o tamanho de cada conjunto e com Jaccard 
# igual a 1.
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
# 5. Guardar como TSV, com o nome do GWAS na primeira coluna
# ----------------------------------------------------------------------------

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
# 6. Estatistica
# ----------------------------------------------------------------------------

cat("\n=== Estatistica ===\n")

# Para a matriz Jaccard (uniao), top 10 pares mais semelhantes
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

cat("\nTop 10 pares de GWAS mais semelhantes (uniao, Jaccard):\n")
print(head(pares, 10))

cat("\nTop 10 pares de GWAS menos semelhantes (uniao, Jaccard maior que zero):\n")
pares_nao_zero <- pares[pares$Jaccard > 0, ]
print(tail(pares_nao_zero, 10))

# Guardar o ranking de pares
write.table(pares, "ranking_pares_gwas_uniao.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\nGuardado: ranking_pares_gwas_uniao.tsv\n")

# Guardar tudo para a fase seguinte no Excel
saveRDS(list(matrizes_uniao      = matrizes_uniao,
             matrizes_intersecao = matrizes_intersecao,
             labels_gwas         = labels_gwas,
             ids_gwas            = ids_gwas_ordenados,
             ranking_pares       = pares),
        "matrizes_gwas.rds")
cat("Guardado: matrizes_gwas.rds (para uso no Excel)\n")