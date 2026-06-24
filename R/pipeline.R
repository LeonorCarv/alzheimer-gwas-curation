# ============================================================================
# pipeline.R
#
# Pipeline de processamento do projeto, desde a extracao das
# plataformas ate a geracao do ficheiro Excel suplementar.
#
# Pre-requisitos (ficheiros que devem estar na pasta antes de correr):
#
#   Manifestos UCSC (descarregados manualmente):
#     - snpArrayAffy250Nsp.txt
#     - snpArrayAffy250Sty.txt
#     - snpArrayAffy6.txt
#     - snpArrayIllumina1M.txt
#     - snpArrayIlluminaHuman660W_Quad.txt
#     - snpArrayIllumina300.txt
#     - snpArrayIllumina550.txt
#     - snpArrayIllumina650.txt
#     - snpArrayIlluminaHumanOmni1_Quad.txt
#
#   Axiom (Thermo Fisher):
#     - Axiom_GW_Hu_SNP.na35.annot.sqlite
#
#   Manifestos Illumina (descarregados manualmente):
#     - GSA-24v1-0_C1.csv
#     - HumanCoreExome-12-v1-0-D.csv + HumanCoreExome-12-v1-0-D-auxilliary-file.txt
#     - HumanCoreExome-12-v1-1-C.csv + HumanCoreExome-12-v1-1-C-auxilliary-file.txt
#     - HumanExome-12-v1-0-B.csv
#     - HumanExome-12-v1-1-B.csv
#     - HumanOmni2.5-4v1_H.csv (humanomni2.5-4v1_h.csv)
#     - HumanOmniExpress-12-v1-0-K.csv
#     - HumanOmniExpress-24-v1-0-B.csv
#     - HumanOmniExpress-24v1-1_A.annotated.txt (Annotation File)
#     - HumanOmniExpressExome-8-v1-2-B.csv
#     - Multi-EthnicGlobal_D1.csv
#     - PsychArray_A_annotated.txt + PsychArray-B-auxiliary-file.txt
#     - Human610-Quadv1_H.bed
#
#   Dados de GWAS:
#     - metadata.csv
#
# Para correr no R, basta fazer:
#   source("pipeline_completo.R")
# ============================================================================

# Limpar sessao
rm(list = ls())

cat("\n###################################################################\n")
cat("# PIPELINE COMPLETO - INICIO\n")
cat("###################################################################\n\n")
cat("ATENCAO: este pipeline e demorado (15-30 minutos).\n")
cat("Para mudancas apenas no metadata.csv, usa pipeline.R em vez deste.\n\n")

inicio_global <- Sys.time()

# ----------------------------------------------------------------------------
# FASE 1: Extracao das plataformas
# Le manifestos e produz padrao_<plataforma>.tsv para cada array
# ----------------------------------------------------------------------------

cat("\n###################################################################\n")
cat("# FASE 1: EXTRACAO DAS PLATAFORMAS\n")
cat("###################################################################\n\n")

source("funcoes_padronizacao.R")

inicio_fase <- Sys.time()

# Bloco UCSC
cat("\n>>> 1.1 Affymetrix (UCSC) <<<\n\n")
processar_ucsc("snpArrayAffy250Nsp.txt.gz",            "Affy250K-Nsp")
processar_ucsc("snpArrayAffy250Sty.txt.gz",            "Affy250K-Sty")
processar_ucsc("snpArrayAffy6.txt.gz",                 "AffySNP6")

cat("\n>>> 1.2 Illumina (UCSC) <<<\n\n")
processar_ucsc("snpArrayIllumina1M.txt.gz",            "Illumina1M")
processar_ucsc("snpArrayIlluminaHuman660W_Quad.txt.gz","Human660W-Quad")
processar_ucsc("snpArrayIllumina300.txt.gz",           "HumanHap300")
processar_ucsc("snpArrayIllumina550.txt.gz",           "HumanHap550")
processar_ucsc("snpArrayIllumina650.txt.gz",           "HumanHap650Y")
processar_ucsc("snpArrayIlluminaHumanOmni1_Quad.txt.gz","HumanOmni1-Quad")

cat("\n>>> 1.3 Thermo Fisher Axiom <<<\n\n")
processar_axiom("Axiom_GW_Hu_SNP.na35.annot.db", "Axiom-CEU")

cat("\n>>> 1.4 Illumina (manifesto CSV padrao) <<<\n\n")
processar_illumina("GSA-24v1-0_C1.csv",              "GSA-v1")
processar_illumina("HumanExome-12-v1-0-B.csv",       "HumanExome-12-v1-0")
processar_illumina("HumanExome-12-v1-1-B.csv",       "HumanExome-12-v1-1")
processar_illumina("HumanOmniExpress-12-v1-0-K.csv", "HumanOmniExpress-12-v1-0")
processar_illumina("HumanOmniExpress-24-v1-0-B.csv", "HumanOmniExpress-24-v1-0")
processar_illumina("HumanOmniExpressExome-8-v1-2-B.csv", "HumanOmniExpressExome-8-v1-2")
processar_illumina("humanomni2.5-4v1_h.csv",         "HumanOmni2.5-Quad")
processar_illumina("Multi-EthnicGlobal_D1.csv",      "MEGA")

cat("\n>>> 1.5 Illumina (manifesto CSV + ficheiro auxiliar) <<<\n\n")
processar_illumina_com_aux(
  ficheiro     = "HumanCoreExome-12-v1-0-D.csv",
  ficheiro_aux = "HumanCoreExome-12-v1-0-D-auxilliary-file.txt",
  plataforma   = "HumanCoreExome-12-v1-0"
)
processar_illumina_com_aux(
  ficheiro     = "HumanCoreExome-12-v1-1-C.csv",
  ficheiro_aux = "HumanCoreExome-12-v1-1-C-auxilliary-file.txt",
  plataforma   = "HumanCoreExome-12-v1-1"
)

cat("\n>>> 1.6 Illumina (Annotation File) <<<\n\n")
processar_illumina_anot(
  ficheiro   = "HumanOmniExpress-24v1-1_A.annotated.txt",
  plataforma = "HumanOmniExpress-24-v1-1"
)

cat("\n>>> 1.7 Illumina (Annotation File + ficheiro auxiliar) <<<\n\n")
processar_illumina_anot_com_aux(
  ficheiro     = "PsychArray_A_annotated.txt",
  ficheiro_aux = "PsychArray-B-auxiliary-file.txt",
  plataforma   = "PsychArray"
)

cat("\n>>> 1.8 Illumina (BED file) <<<\n\n")
processar_bed("Human610-Quadv1_H.bed", "Human610-Quad")

fim_fase <- Sys.time()
cat("\nDuracao da FASE 1:",
    round(as.numeric(difftime(fim_fase, inicio_fase, units = "mins")), 2),
    "minutos\n")

# ----------------------------------------------------------------------------
# FASE 2: Juncao das plataformas num ficheiro mestre
# Combina todos os padrao_*.tsv num unico padrao_TODOS.tsv
# ----------------------------------------------------------------------------

cat("\n###################################################################\n")
cat("# FASE 2: JUNCAO DAS PLATAFORMAS\n")
cat("###################################################################\n\n")

inicio_fase <- Sys.time()

mestre <- juntar_padroes()

fim_fase <- Sys.time()
cat("\nDuracao da FASE 2:",
    round(as.numeric(difftime(fim_fase, inicio_fase, units = "mins")), 2),
    "minutos\n")

# ----------------------------------------------------------------------------
# FASE 3: Analise de cruzamento entre plataformas
# Matriz e distribuicao de SNPs partilhados entre plataformas
# ----------------------------------------------------------------------------

cat("\n###################################################################\n")
cat("# FASE 3: ANALISE DE CRUZAMENTO ENTRE PLATAFORMAS\n")
cat("###################################################################\n\n")

inicio_fase <- Sys.time()

source("analise_cruzamento_platf.R")

fim_fase <- Sys.time()
cat("\nDuracao da FASE 3:",
    round(as.numeric(difftime(fim_fase, inicio_fase, units = "mins")), 2),
    "minutos\n")

# ----------------------------------------------------------------------------
# FASE 4: Analise por GWAS
# Para cada GWAS, calcular SNPs por uniao e intersecao
# ----------------------------------------------------------------------------

cat("\n###################################################################\n")
cat("# FASE 4: ANALISE POR GWAS\n")
cat("###################################################################\n\n")

inicio_fase <- Sys.time()

source("preparar_gwas.R")
source("analise_gwas.R")
source("analise_gwas_cruzamento.R")

fim_fase <- Sys.time()
cat("\nDuracao da FASE 4:",
    round(as.numeric(difftime(fim_fase, inicio_fase, units = "mins")), 2),
    "minutos\n")

# ----------------------------------------------------------------------------
# FASE 5: Geracao do ficheiro Excel suplementar
# ----------------------------------------------------------------------------

cat("\n###################################################################\n")
cat("# FASE 5: FICHEIRO EXCEL SUPLEMENTAR\n")
cat("###################################################################\n\n")

inicio_fase <- Sys.time()

source("generate_supplementary_excel.R")

fim_fase <- Sys.time()
cat("\nDuracao da FASE 5:",
    round(as.numeric(difftime(fim_fase, inicio_fase, units = "mins")), 2),
    "minutos\n")

# ----------------------------------------------------------------------------
# Fim
# ----------------------------------------------------------------------------
fim_global <- Sys.time()
duracao <- difftime(fim_global, inicio_global, units = "mins")

cat("\n###################################################################\n")
cat("# PIPELINE COMPLETO - FIM\n")
cat("###################################################################\n\n")
cat("Duracao total:", round(as.numeric(duracao), 2), "minutos\n")
cat("Inicio:", format(inicio_global, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Fim:   ", format(fim_global,    "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Ficheiros chave gerados:\n")
ficheiros_chave <- c(
  # Padroes por plataforma
  "padrao_Affy250K-Nsp.tsv",
  "padrao_Affy250K-Sty.tsv",
  "padrao_AffySNP6.tsv",
  "padrao_Axiom-CEU.tsv",
  "padrao_GSA-v1.tsv",
  "padrao_Human610-Quad.tsv",
  "padrao_Human660W-Quad.tsv",
  "padrao_HumanCoreExome-12-v1-0.tsv",
  "padrao_HumanCoreExome-12-v1-1.tsv",
  "padrao_HumanExome-12-v1-0.tsv",
  "padrao_HumanExome-12-v1-1.tsv",
  "padrao_HumanHap300.tsv",
  "padrao_HumanHap550.tsv",
  "padrao_HumanHap650Y.tsv",
  "padrao_HumanOmni1-Quad.tsv",
  "padrao_HumanOmni2.5-Quad.tsv",
  "padrao_HumanOmniExpress-12-v1-0.tsv",
  "padrao_HumanOmniExpress-24-v1-0.tsv",
  "padrao_HumanOmniExpress-24-v1-1.tsv",
  "padrao_HumanOmniExpressExome-8-v1-2.tsv",
  "padrao_Illumina1M.tsv",
  "padrao_MEGA.tsv",
  "padrao_PsychArray.tsv",
  # Mestre e analise de cruzamento entre plataformas
  "padrao_TODOS.tsv",
  "resumo_plataformas.tsv",
  "matriz_sobreposicao.tsv",
  "snps_por_n_plataformas.tsv",
  # Analise por GWAS
  "gwas_limpo.tsv",
  "gwas_snps_sumario.tsv",
  "conjuntos_snps_por_gwas.rds",
  "matriz_gwas_uniao_absoluta.tsv",
  "matriz_gwas_uniao_jaccard.tsv",
  "matriz_gwas_intersecao_absoluta.tsv",
  "matriz_gwas_intersecao_jaccard.tsv",
  "ranking_pares_gwas_uniao.tsv",
  "matrizes_gwas.rds",
  # Excel final
  "Supplementary_tables_genotyping_platforms.xlsx"
)
n_ok <- 0
n_falha <- 0
for (f in ficheiros_chave) {
  if (file.exists(f)) {
    n_ok <- n_ok + 1
    cat(sprintf("  [OK]  %s\n", f))
  } else {
    n_falha <- n_falha + 1
    cat(sprintf("  [!!]  %s (nao encontrado)\n", f))
  }
}
cat(sprintf("\nResumo: %d ficheiros gerados, %d em falta.\n\n",
            n_ok, n_falha))