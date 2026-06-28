# ============================================================================
# preparar_gwas.R
# Estrutura a tabela de GWAS para a analise de cruzamento.
#
# Regras importantes no tratamento das celulas:
#   - Célula vazia no CSV (sem conteúdo) era uma célula merged no Excel e herda
#     o valor da célula imediatamente acima, dentro do mesmo GWAS
#   - Célula com "NA" textual é uma plataforma desconhecida no paper original, 
#     e mantém-se como NA.
# ============================================================================

library(zoo)

# ----------------------------------------------------------------------------
# 1. Ler o CSV
# na.strings = "" apenas, para distinguir celulas vazias (merged)
# das células com "NA" textual (desconhecidas)
# ----------------------------------------------------------------------------

gwas <- read.csv2("../data/raw/metadata.csv",
                  stringsAsFactors = FALSE,
                  na.strings = "",
                  check.names = FALSE,
                  fileEncoding = "Windows-1252")

cat("Total de linhas lidas:", nrow(gwas), "\n")
cat("Total de colunas lidas:", ncol(gwas), "\n\n")

# Remover linhas inteiramente vazias (linhas em branco no fim do CSV)
linhas_vazias <- is.na(gwas$`Studies where the data comes from`)
n_vazias <- sum(linhas_vazias)
gwas <- gwas[!linhas_vazias, ]
cat("Linhas vazias removidas:", n_vazias, "\n")
cat("Linhas restantes apos remover vazias:", nrow(gwas), "\n\n")

# ----------------------------------------------------------------------------
# 2. Preenchimento descendente dos metadados do GWAS
# Estes campos só estão preenchidos na primeira linha de cada GWAS
# ----------------------------------------------------------------------------

gwas$`FIRST AUTHOR` <- na.locf(gwas$`FIRST AUTHOR`, na.rm = FALSE)
gwas$`Year of publication` <- na.locf(gwas$`Year of publication`, na.rm = FALSE)
gwas$ID <- na.locf(gwas$ID, na.rm = FALSE)

if ("Web/ citation" %in% colnames(gwas)) {
  gwas$`Web/ citation` <- na.locf(gwas$`Web/ citation`, na.rm = FALSE)
}
if ("Title" %in% colnames(gwas)) {
  gwas$Title <- na.locf(gwas$Title, na.rm = FALSE)
}

# ----------------------------------------------------------------------------
# 3. Preenchimento descendente da coluna Genotyping Platform dentro de cada 
# GWAS, mantendo "NA" textual como genuinamente desconhecida
#   - Célula vazia (NA no R) herda da célula acima
#   - Célula com "NA" textual não herda e fica como string "NA"
# Após o preenchimento, o "NA" textual passa a NA real do R
# ----------------------------------------------------------------------------

cat("Antes do preenchimento descendente:\n")
cat("  Células vazias na plataforma:", sum(is.na(gwas$`Genotyping Platform`)), "\n")
cat("  Células com 'NA' textual:", sum(gwas$`Genotyping Platform` == "NA",
                                       na.rm = TRUE), "\n\n")

# Aplicar na.locf por GWAS, para nao herdar plataformas de um GWAS para outro
gwas_lista <- split(gwas, gwas$ID)
gwas_lista <- lapply(gwas_lista, function(grupo) {
  grupo$`Genotyping Platform` <- na.locf(grupo$`Genotyping Platform`,
                                         na.rm = FALSE)
  grupo$`Genomic Assay GWAS` <- na.locf(grupo$`Genomic Assay GWAS`,
                                        na.rm = FALSE)
  grupo
})
gwas <- do.call(rbind, gwas_lista)
rownames(gwas) <- NULL

# Converter "NA" textual para NA real do R (subamostras genuinamente desconhecidas)
gwas$`Genotyping Platform`[gwas$`Genotyping Platform` == "NA"] <- NA

cat("Após preenchimento descendente:\n")
cat("  Células com plataforma definida:",
    sum(!is.na(gwas$`Genotyping Platform`)), "\n")
cat("  Células com plataforma NA (genuinamente desconhecida):",
    sum(is.na(gwas$`Genotyping Platform`)), "\n\n")

# ----------------------------------------------------------------------------
# 3. Remover linhas TOTAL
# ----------------------------------------------------------------------------
n_antes <- nrow(gwas)
gwas <- gwas[gwas$`Studies where the data comes from` != "TOTAL", ]
cat("Linhas TOTAL removidas:", n_antes - nrow(gwas), "\n")
cat("Linhas restantes (subamostras):", nrow(gwas), "\n\n")

# ----------------------------------------------------------------------------
# 4. Renomear colunas para nomes simples
# ----------------------------------------------------------------------------
colnames(gwas)[colnames(gwas) == "FIRST AUTHOR"] <- "autor"
colnames(gwas)[colnames(gwas) == "Year of publication"] <- "ano"
colnames(gwas)[colnames(gwas) == "Studies where the data comes from"] <- "subamostra"
colnames(gwas)[colnames(gwas) == "N - AD"] <- "n_ad"
colnames(gwas)[colnames(gwas) == "N - Controls"] <- "n_controlos"
colnames(gwas)[colnames(gwas) == "N - Total"] <- "n_total"
colnames(gwas)[colnames(gwas) == "Genotyping Platform"] <- "plataformas"
colnames(gwas)[colnames(gwas) == "Genomic Assay GWAS"] <- "genomic_assay"

# ----------------------------------------------------------------------------
# 5. Diagnostico final
# ----------------------------------------------------------------------------
cat("GWAS distintos:", length(unique(gwas$ID)), "\n")

cat("\nSubamostras por GWAS:\n")
print(table(gwas$ID))

cat("\nSubamostras com plataforma definida (válidas):",
    sum(!is.na(gwas$plataformas)), "\n")
cat("Subamostras com plataforma NA (desconhecidas):",
    sum(is.na(gwas$plataformas)), "\n\n")

# Distribuição do número de plataformas por subamostra
n_plats <- sapply(strsplit(gwas$plataformas, ";"), length)
cat("Distribuicao do numero de plataformas por subamostra:\n")
print(table(n_plats, useNA = "ifany"))

# ----------------------------------------------------------------------------
# 6. Mostrar primeiras linhas
# ----------------------------------------------------------------------------
cat("\nPrimeiras 20 linhas apos limpeza:\n")
print(head(gwas[, c("autor", "ano", "ID", "subamostra",
                    "n_ad", "n_controlos", "plataformas")], 20))

# ----------------------------------------------------------------------------
# 7. Guardar versão limpa
# ----------------------------------------------------------------------------
gwas_essencial <- gwas[, c("autor", "ano", "ID", "subamostra",
                           "n_ad", "n_controlos", "n_total", "plataformas",
                           "genomic_assay")]
write.table(gwas_essencial, "gwas_limpo.tsv", sep = "\t",
            row.names = FALSE, quote = TRUE, na = "NA")
cat("\nGuardado: gwas_limpo.tsv\n")
