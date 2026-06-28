# ============================================================================
# analise_cruzamento_platf.R
#
# Analise de cruzamento entre plataformas a partir do ficheiro mestre
# padrao_TODOS.tsv.
#
# A chave de cruzamento é (chr, pos), sem o identificador rs, para evitar
# problemas de harmonizacao entre versoes do dbSNP (UCSC build 131-132,
# Axiom build 142, Illumina variavel).
#
# Saidas:
#   - matriz_sobreposicao.tsv: matriz quadrada com SNPs partilhados
#     entre cada par de plataformas
#   - resumo_plataformas.tsv: por plataforma, total de SNPs unicos
#     (apos remocao de duplicados na propria plataforma)
#   - snps_por_n_plataformas.tsv: distribuicao do numero de plataformas
#     que cobrem cada SNP unico
# ============================================================================

mestre <- read.table("padrao_TODOS.tsv",
                     sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE,
                     na.strings = "NA", quote = "", comment.char = "")

cat("Total de linhas no ficheiro mestre:", nrow(mestre), "\n")
cat("Plataformas distintas:", length(unique(mestre$plataforma)), "\n\n")

# ----------------------------------------------------------------------------
# 1. Construir chave (chr, pos) e remover duplicados por plataforma
# (mesma posição registada mais que uma vez na mesma plataforma)
# ----------------------------------------------------------------------------
mestre$chave <- paste(mestre$chr, mestre$pos, sep = ":")

n_antes <- nrow(mestre)
mestre <- mestre[!duplicated(mestre[, c("plataforma", "chave")]), ]
n_depois <- nrow(mestre)
cat("Duplicados intra-plataforma removidos:", n_antes - n_depois, "\n\n")

# ----------------------------------------------------------------------------
# 2. Resumo por plataforma: Nº de SNPs unicos por plataforma
# ----------------------------------------------------------------------------
resumo <- aggregate(chave ~ plataforma, data = mestre, FUN = length)
colnames(resumo) <- c("plataforma", "n_snps_unicos")
resumo <- resumo[order(-resumo$n_snps_unicos), ]
print(resumo)

write.table(resumo, "resumo_plataformas.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------------------------------------------------------
# 3. Matriz de sobreposicao entre plataformas
# ----------------------------------------------------------------------------

# Construir lista de SNPs (chaves) por plataforma
plataformas <- sort(unique(mestre$plataforma))
n_plat <- length(plataformas)
snps_por_plat <- split(mestre$chave, mestre$plataforma)

# Matriz quadrada n_plat x n_plat
matriz <- matrix(0L, nrow = n_plat, ncol = n_plat,
                 dimnames = list(plataformas, plataformas))

for (i in seq_along(plataformas)) {
  snps_i <- snps_por_plat[[plataformas[i]]]
  for (j in seq_along(plataformas)) {
    if (i == j) {
      matriz[i, j] <- length(snps_i)
    } else if (i < j) {
      snps_j <- snps_por_plat[[plataformas[j]]]
      n_partilhados <- length(intersect(snps_i, snps_j))
      matriz[i, j] <- n_partilhados
      matriz[j, i] <- n_partilhados
    }
  }
  cat(sprintf("  %d / %d - %s\n", i, n_plat, plataformas[i]))
}

write.table(matriz, "matriz_sobreposicao.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# ----------------------------------------------------------------------------
# 4. SNPs partilhados entre N plataformas
# ----------------------------------------------------------------------------

# Para cada SNP unico (chave), contar em quantas plataformas aparece
contagem_por_snp <- aggregate(plataforma ~ chave, data = mestre,
                              FUN = function(x) length(unique(x)))
colnames(contagem_por_snp) <- c("chave", "n_plataformas")

distribuicao <- table(contagem_por_snp$n_plataformas)
distribuicao_df <- data.frame(
  n_plataformas = as.integer(names(distribuicao)),
  n_snps        = as.integer(distribuicao)
)
distribuicao_df <- distribuicao_df[order(distribuicao_df$n_plataformas), ]
print(distribuicao_df)

write.table(distribuicao_df, "snps_por_n_plataformas.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Total de SNPs únicos (após união de todas as plataformas):",
    nrow(contagem_por_snp), "\n")
cat("SNPs presentes em apenas uma plataforma:",
    sum(contagem_por_snp$n_plataformas == 1), "\n")
cat("SNPs presentes em duas ou mais plataformas:",
    sum(contagem_por_snp$n_plataformas >= 2), "\n")
cat("SNPs presentes em todas as plataformas:",
    sum(contagem_por_snp$n_plataformas == n_plat), "\n")