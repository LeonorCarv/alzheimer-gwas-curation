# ============================================================================
# analise_gwas.R
#
# Calcula, para cada GWAS:
#   - n_snps_uniao: numero de SNPs na UNIAO de todas as subamostras com
#     plataforma conhecida. Conjunto otimista (cobertura potencial).
#   - n_snps_intersecao: numero de SNPs na INTERSECAO das unioes de cada
#     subamostra com plataforma conhecida. Conjunto biologicamente analisavel
#     antes de imputacao.
#
# Regras de agregação:
#   - Subamostra com plataforma NA ou nao mapeada:  ignorada
#   - Subamostra com multiplas plataformas:         UNIAO dos SNPs das plataformas
#   - GWAS com multiplas subamostras:               INTERSECAO das unioes de cada subamostra
#
# Lê gwas_limpo.tsv e padrao_TODOS.tsv (caminhos relativos a raiz do projeto).
# ============================================================================

# ----------------------------------------------------------------------------
# 1. Ler dados
# ----------------------------------------------------------------------------

gwas <- read.table("gwas_limpo.tsv", sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE, na.strings = "NA",
                   quote = "\"", comment.char = "")

mestre <- read.table("padrao_TODOS.tsv", sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE, na.strings = "NA",
                     quote = "", comment.char = "")

mestre$chave <- paste(mestre$chr, mestre$pos, sep = ":")
mestre <- mestre[!duplicated(mestre[, c("plataforma", "chave")]), ]

cat("GWAS distintos:", length(unique(gwas$ID)), "\n")
cat("Subamostras no total:", nrow(gwas), "\n")

# ----------------------------------------------------------------------------
# 2. Mapeamento do nome da plataforma no CSV para o nome do ficheiro padrao
# Agrupado por origem do manifesto.
# ----------------------------------------------------------------------------

mapeamento <- c(
  # Affymetrix UCSC
  "Affymetrix GeneChip Human Mapping 250K Nsp Array"    = "Affy250K-Nsp",
  "Affymetrix GeneChip Human Mapping 250K Sty Array"    = "Affy250K-Sty",
  "Applied Biosystems Human Genome-Wide SNP Array 6.0"  = "AffySNP6",
  
  # Thermo Fisher Axiom
  "Applied Biosystems Axiom Genome-Wide CEU Array"      = "Axiom-CEU",
  
  # Illumina UCSC
  "Illumina Infinium Human1M-Duo"                       = "Illumina1M",
  "Illumina Infinium Human660W-Quad"                    = "Human660W-Quad",
  "Illumina Infinium HumanHap300"                       = "HumanHap300",
  "Illumina Infinium HumanHap550"                       = "HumanHap550",
  "Illumina Infinium HumanHap650Y"                      = "HumanHap650Y",
  "Illumina Infinium HumanOmni1-Quad"                   = "HumanOmni1-Quad",
  
  # Illumina manifesto direto
  "Illumina Infinium Global Screening Array"            = "GSA-v1",
  "Illumina Infinium Human610-Quad"                     = "Human610-Quad",
  "Illumina Infinium HumanCoreExome-12 v1.0"            = "HumanCoreExome-12-v1-0",
  "Illumina Infinium HumanCoreExome12 v1.1"             = "HumanCoreExome-12-v1-1",
  "Illumina Infinium HumanExome v1.0"                   = "HumanExome-12-v1-0",
  "Illumina Infinium HumanExome v1.1"                   = "HumanExome-12-v1-1",
  "Illumina Infinium HumanOmni2.5-Quad"                 = "HumanOmni2.5-Quad",
  "Illumina Infinium HumanOmniExpress-12"               = "HumanOmniExpress-12-v1-0",
  "Illumina Infinium HumanOmniExpress-24 v1.0"          = "HumanOmniExpress-24-v1-0",
  "Illumina Infinium OmniExpress-24 v1.0"               = "HumanOmniExpress-24-v1-0",
  "Illumina Infinium HumanOmniExpress-24 v1.1"          = "HumanOmniExpress-24-v1-1",
  "Illumina Infinium HumanOmniExpressExome-8 v1.2"      = "HumanOmniExpressExome-8-v1-2",
  "Illumina Infinium Multi-Ethnic Global-8 Array (MEGA)"= "MEGA",
  "Illumina Infinium Expanded Multi-Ethnic Genotyping Array (MEGAEX)" = "MEGA",
  "Illumina Infinium PsychArray"                        = "PsychArray"
)

# Plataformas excluídas (sem ficheiro padrao - as subamostras que so as usam são ignoradas)
plataformas_excluidas <- c(
  "Agena Bioscience MassARRAY System",
  "Applied Biosystems Axiom Spain Biobank Array",
  "Applied Biosystems UK BiLEVE Axiom Array",
  "Applied Biosystems UK Biobank Axiom Array",
  "deCODE genetics Custom Axiom Array v1 (Part No. 20012591_A1)",
  "Illumina GoldenGate Custom Genotyping",
  "Illumina HiSeq X Ten",
  "Illumina Infinium HumanCNV370-Duo",
  "Illumina Infinium Metabochip",
  "Illumina Infinium Neurochip",
  "Illumina Infinium UM HUNT Biobank v1.0",
  "Illumina iSelect Custom Genotyping BeadChip",
  "Illumina NovaSeq 6000",
  "LGC KASP Genotyping Assay",
  "Sequenom MassARRAY System",
  "TaqMan Genotyping Assay"
)

# ----------------------------------------------------------------------------
# 3. Funcao auxiliar: Devolve o nome do ficheiro padrao para um nome do CSV, 
# ou NA se a plataforma for excluída ou não estiver mapeada
# ----------------------------------------------------------------------------

mapear_plataforma <- function(nome_csv) {
  nome_csv <- trimws(nome_csv)
  if (is.na(nome_csv) || nome_csv == "") return(NA_character_)
  if (nome_csv %in% plataformas_excluidas) return(NA_character_)
  if (nome_csv %in% names(mapeamento)) return(mapeamento[[nome_csv]])
  return(NA_character_)  # nao reconhecida
}

# ----------------------------------------------------------------------------
# 4. Construir lista de SNPs por plataforma (chaves chr:pos)
# ----------------------------------------------------------------------------

snps_por_plat <- split(mestre$chave, mestre$plataforma)

# ----------------------------------------------------------------------------
# 5. Calcular conjunto de SNPs por subamostra (união das plataformas)
# ----------------------------------------------------------------------------

gwas$plataformas_mapeadas <- ""
gwas$snps_uniao_subamostra <- I(vector("list", nrow(gwas)))
gwas$n_plataformas_validas <- 0L
gwas$subamostra_valida <- FALSE

for (i in seq_len(nrow(gwas))) {
  plats_csv <- gwas$plataformas[i]
  if (is.na(plats_csv) || plats_csv == "") {
    gwas$snps_uniao_subamostra[[i]] <- character(0)
    next
  }
  
  # Separar plataformas (por ponto-e-virgula)
  plats_lista <- trimws(strsplit(plats_csv, ";")[[1]])
  
  # Mapear cada uma
  plats_mapeadas <- sapply(plats_lista, mapear_plataforma, USE.NAMES = FALSE)
  plats_validas <- plats_mapeadas[!is.na(plats_mapeadas)]
  
  gwas$plataformas_mapeadas[i] <- paste(plats_validas, collapse = "; ")
  gwas$n_plataformas_validas[i] <- length(plats_validas)
  
  if (length(plats_validas) == 0) {
    gwas$snps_uniao_subamostra[[i]] <- character(0)
    next
  }
  
  # União dos SNPs das plataformas válidas desta subamostra
  snps_lista <- snps_por_plat[plats_validas]
  snps_uniao <- unique(unlist(snps_lista))
  
  gwas$snps_uniao_subamostra[[i]] <- snps_uniao
  gwas$subamostra_valida[i] <- TRUE
}

cat("Subamostras com plataformas válidas:", sum(gwas$subamostra_valida), "\n")
cat("Subamostras ignoradas (NA ou excluidas):", sum(!gwas$subamostra_valida), "\n")

# ----------------------------------------------------------------------------
# Calcular conjunto de SNPs por GWAS (agregar por GWAS)
#   - União junta todos os SNPs de qualquer subamostra válida
#   - Interseção guarda os SNPs presentes em todas elas
# ----------------------------------------------------------------------------

ids_gwas <- unique(gwas$ID)
resultado <- data.frame(
  ID                            = character(),
  autor                         = character(),
  ano                           = integer(),
  n_subamostras_total           = integer(),
  n_subamostras_com_plataforma  = integer(),
  n_subamostras_NA              = integer(),
  plataformas_unicas            = character(),
  n_snps_uniao                  = integer(),
  n_snps_intersecao             = integer(),
  stringsAsFactors              = FALSE
)

for (id in ids_gwas) {
  linhas_gwas <- gwas[gwas$ID == id, ]
  
  n_total <- nrow(linhas_gwas)
  n_validas <- sum(linhas_gwas$subamostra_valida)
  n_na <- n_total - n_validas
  
  # Plataformas únicas usadas neste GWAS (válidas)
  plats_validas_gwas <- unique(unlist(strsplit(
    linhas_gwas$plataformas_mapeadas[linhas_gwas$subamostra_valida],
    "; ")))
  plats_validas_gwas <- sort(plats_validas_gwas)
  
  # União: todos os SNPs de qualquer subamostra válida
  if (n_validas == 0) {
    snps_uniao <- character(0)
    snps_intersecao <- character(0)
  } else {
    snps_uniao <- unique(unlist(linhas_gwas$snps_uniao_subamostra[
      linhas_gwas$subamostra_valida]))
    
    # Interseção: SNPs presentes em todas as subamostras válidas
    snps_intersecao <- Reduce(intersect,
                              linhas_gwas$snps_uniao_subamostra[
                                linhas_gwas$subamostra_valida])
  }
  
  resultado <- rbind(resultado, data.frame(
    ID                            = id,
    autor                         = linhas_gwas$autor[1],
    ano                           = linhas_gwas$ano[1],
    n_subamostras_total           = n_total,
    n_subamostras_com_plataforma  = n_validas,
    n_subamostras_NA              = n_na,
    plataformas_unicas            = paste(plats_validas_gwas, collapse = "; "),
    n_snps_uniao                  = length(snps_uniao),
    n_snps_intersecao             = length(snps_intersecao),
    stringsAsFactors              = FALSE
  ))
}

# Ordenar por ID numérico
resultado$ID_num <- as.integer(sub("GWAS", "", resultado$ID))
resultado <- resultado[order(resultado$ID_num), ]
resultado$ID_num <- NULL

# ----------------------------------------------------------------------------
# 7. Apresentar e gravar resultado
# ----------------------------------------------------------------------------

print(resultado)

write.table(resultado, "gwas_snps_sumario.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE, na = "NA")
cat("\nGuardado: gwas_snps_sumario.tsv\n")

# Estatistica final
cat("\n=== Estatística final ===\n")
cat("GWAS analisados:", nrow(resultado), "\n")
cat("GWAS sem subamostras válidas:",
    sum(resultado$n_subamostras_com_plataforma == 0), "\n")
cat("GWAS com 1+ subamostras válidas:",
    sum(resultado$n_subamostras_com_plataforma >= 1), "\n")

# Guardar os conjuntos de SNPs por GWAS para a fase de cruzamento
conjuntos_snps <- list()
for (id in resultado$ID) {
  linhas_gwas <- gwas[gwas$ID == id & gwas$subamostra_valida, ]
  if (nrow(linhas_gwas) == 0) {
    conjuntos_snps[[id]] <- list(uniao = character(0),
                                 intersecao = character(0))
    next
  }
  uniao <- unique(unlist(linhas_gwas$snps_uniao_subamostra))
  intersecao <- Reduce(intersect, linhas_gwas$snps_uniao_subamostra)
  conjuntos_snps[[id]] <- list(uniao = uniao, intersecao = intersecao)
}

saveRDS(conjuntos_snps, "conjuntos_snps_por_gwas.rds")
cat("Guardado: conjuntos_snps_por_gwas.rds (para uso na fase de cruzamento)\n")