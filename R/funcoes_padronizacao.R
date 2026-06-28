# ============================================================================
# funcoes_padronizacao.R
# Funcoes para padronizar manifestos e tabelas de plataformas de genotipagem
#
# Versao: 1.13
# Ultima alteracao: adicionada processar_bed para ficheiros BED de UCSC/Illumina
#                   (4 colunas: chr, start, end, name), com conversao chr<N> para
#                   <N>, M para MT e zero-based para one-based. O BED nao traz
#                   coluna de alelos, por isso os indels nao sao filtrados neste
#                   caso. Usada para o Human610-Quad apos falha do FTP do manifesto.
#
# Historico:
#   1.0 - Versao inicial com processar_ucsc, processar_illumina, processar_axiom
#   1.1 - Detecao automatica Illumina/Affymetrix em processar_ucsc
#   1.2 - criar_padrao passa a aceitar marcadores posicionais (GA, VG, etc.)
#         excluindo apenas cnvi e Mito por prefixo
#   1.3 - criar_padrao passa a filtrar indels e variantes nao-SNP pela coluna
#         observed (aceita apenas formato A/C/G/T)
#   1.4 - processar_ucsc passa a converter posicoes de zero-based para
#         one-based, em coerencia com a convencao usada pela Illumina e
#         pela Thermo Fisher (Axiom)
#   1.5 - processar_axiom adaptada a estrutura real dos ficheiros Axiom:
#         conversao de cromossomas Affymetrix (24=X, 25=Y, 26=MT), construcao
#         de pseudo-coluna observed a partir de Allele_A e Allele_B para
#         permitir o filtro de indels, e validacao contra ficheiros fantasma
#   1.6 - processar_illumina adaptada para o formato real dos manifestos
#         Illumina: corte de [Controls], conversao do cromossoma XY para X
#         (regiao pseudoautosomica), e passagem da coluna SNP a criar_padrao
#         no formato A/B para filtrar indels (I/D, D/I)
#   1.8 - adicionada funcao processar_illumina_com_aux para integrar ficheiros
#         auxiliares de mapeamento Name -> RsID em arrays Illumina especificos
#         (CoreExome v1.0 e v1.1). Aplica recuperacao de rs canonico apenas
#         para Names nao-rs mapeados para um unico rs valido. Casos com rs
#         antigos para rs novos, multiplos rs, ou rs descontinuados sao
#         rejeitados em coerencia com a regra de deferir harmonizacao
#         entre versoes do dbSNP para passo posterior.
#   1.9 - adicionada funcao processar_illumina_anot para processar
#         Annotation Files da Illumina (formato Name, Chr, MapInfo, Alleles)
#         para arrays sem manifesto CSV padrao acessivel. Assume hg19 como
#         build de referencia. Usada para HumanOmniExpress-24 v1.1 e
#         PsychArray.
#   1.10 - processar_illumina_anot passa a usar fill = TRUE na leitura,
#          para tolerar ficheiros TXT com numero variavel de tabs nas
#          ultimas linhas (caso do HumanOmniExpress-24 v1.1)
#   1.11 - adicionada funcao processar_illumina_anot_com_aux, que combina
#          a leitura de Annotation Files (com fill = TRUE) com a integracao
#          de ficheiros auxiliares de mapeamento Name -> RsID. Usada para
#          o PsychArray, onde cerca de 49 por cento dos marcadores tem
#          mapeamento aplicavel no auxiliar.
#   1.12 - processar_illumina passa a normalizar GenomeBuild "37.x" para
#          "37" antes de criar o ficheiro padrao, em coerencia com os
#          outros arrays anotados como hg19. Necessario para o
#          HumanOmni2.5-Quad, que vem anotado contra "37.1" (GRCh37 com
#          correcoes menores subsequentes).
#   1.13 - adicionada funcao processar_bed para processar ficheiros BED
#          de UCSC/Illumina (4 colunas: chr, start, end, name). Faz
#          conversao chr<N> para <N>, M para MT, e zero-based para
#          one-based (start + 1). Limitacao metodologica: BED nao tem
#          coluna de alelos, indels (se existirem) NAO sao filtrados.
#          Usada para o Human610-Quad apos falha do FTP do manifesto.
# ============================================================================


# ---------------------------------------------------------------------------
# 1. CONFIGURACAO INICIAL
# ---------------------------------------------------------------------------
library(DBI)
library(RSQLite)


# ---------------------------------------------------------------------------
# 2. FUNCAO PRINCIPAL DE PADRONIZACAO
# ---------------------------------------------------------------------------
# Recebe dados em bruto de qualquer fonte e devolve um data frame no formato
# padrao do projeto: plataforma, rs, chr, pos, build, tipo_id, nome_original.
#
# Aplica em conjunto:
#   - Limpeza do prefixo "chr" do cromossoma
#   - Conversao da build numerica (36, 37, 38) para nomenclatura hg
#   - Filtragem de marcadores nao-SNP por prefixo: cnvi (CNV) e Mito/MT
#   - Filtragem de variantes nao-SNP pela coluna observed: aceita apenas A/C/G/T
#     (exclui indels I/D, D/I e outros formatos)
#   - Exclusao de cromossomas nao nucleares (1-22, X, Y apenas)
#   - Classificacao por tipo de identificador (rs ou posicao)
# ---------------------------------------------------------------------------
criar_padrao <- function(rs_bruto, chr_bruto, pos_bruto, build_bruto,
                         plataforma, observed_bruto = NULL) {
  
  rs_bruto    <- as.character(rs_bruto)
  chr_bruto   <- as.character(chr_bruto)
  build_bruto <- as.character(build_bruto)
  pos_bruto   <- as.integer(pos_bruto)
  
  chr_limpo <- sub("^chr", "", chr_bruto)
  
  build_limpa <- build_bruto
  build_limpa[build_limpa == "36"] <- "hg18"
  build_limpa[build_limpa == "37"] <- "hg19"
  build_limpa[build_limpa == "38"] <- "hg38"
  
  # Categorias a EXCLUIR por nao serem SNPs nucleares (por prefixo)
  e_cnvi <- grepl("^cnvi", rs_bruto)
  e_mito <- grepl("^Mito|^MT", rs_bruto)
  excluir_biologico <- e_cnvi | e_mito
  
  # Coordenadas validas
  cromossoma_valido <- chr_limpo %in% c(as.character(1:22), "X", "Y")
  posicao_valida    <- !is.na(pos_bruto) & pos_bruto > 0
  coordenadas_ok    <- cromossoma_valido & posicao_valida
  
  # Filtro de SNP biologico (so se a coluna observed estiver disponivel):
  # aceitar apenas variantes A/C/G/T, excluir indels (I/D, D/I) e outros formatos
  if (!is.null(observed_bruto)) {
    observed_bruto <- as.character(observed_bruto)
    e_snp_real <- grepl("^[ACGT]/[ACGT]$", observed_bruto)
  } else {
    e_snp_real <- rep(TRUE, length(rs_bruto))
  }
  
  e_rs <- grepl("^rs[0-9]+", rs_bruto) & coordenadas_ok &
    !excluir_biologico & e_snp_real
  e_posicional <- coordenadas_ok & !excluir_biologico & !e_rs & e_snp_real
  
  tipo <- ifelse(e_rs, "rs",
                 ifelse(e_posicional, "posicao", "excluir"))
  
  padrao <- data.frame(
    plataforma    = plataforma,
    rs            = ifelse(tipo == "rs", rs_bruto, NA),
    chr           = chr_limpo,
    pos           = pos_bruto,
    build         = build_limpa,
    tipo_id       = tipo,
    nome_original = rs_bruto,
    stringsAsFactors = FALSE
  )
  
  cat("Plataforma:", plataforma, "\n")
  cat("  Total de linhas no ficheiro:    ", length(rs_bruto), "\n")
  cat("  Excluidas (cnvi):               ", sum(e_cnvi), "\n")
  cat("  Excluidas (Mito):               ", sum(e_mito), "\n")
  if (!is.null(observed_bruto)) {
    cat("  Excluidas (nao-SNP, ex: indels): ", sum(!e_snp_real), "\n")
  }
  cat("  Excluidas (coordenadas invalidas):",
      sum(!coordenadas_ok & !excluir_biologico & e_snp_real), "\n")
  
  padrao <- padrao[padrao$tipo_id != "excluir", ]
  
  cat("  Incluidas (rs):                 ", sum(padrao$tipo_id == "rs"), "\n")
  cat("  Incluidas (posicao):            ", sum(padrao$tipo_id == "posicao"), "\n")
  cat("  Total incluido:                 ", nrow(padrao), "\n\n")
  
  return(padrao)
}

# ---------------------------------------------------------------------------
# 3. FUNCAO DE GRAVACAO
# ---------------------------------------------------------------------------
gravar_padrao <- function(padrao, plataforma) {
  nome_ficheiro <- paste0("padrao_", plataforma, ".tsv")
  write.table(padrao,
              file = nome_ficheiro,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE,
              na = "NA",
              fileEncoding = "UTF-8")
  cat("Gravado:", nome_ficheiro, "\n\n")
}

# ===========================================================================
# 4. PROCESSAMENTO POR FONTE
# ===========================================================================
# Cada bloco trata uma fonte: UCSC, manifesto Illumina, anotacao Axiom SQLite.
# A funcao criar_padrao e gravar_padrao sao sempre as mesmas. So muda a leitura.
# ===========================================================================

# ---------------------------------------------------------------------------
# 4.1. ARRAYS DO UCSC (tabelas snpArray*.txt.gz)
# ---------------------------------------------------------------------------
# As tabelas do UCSC tem duas estruturas conforme o fabricante do array.
# A funcao processar_ucsc deteta automaticamente o tipo pelo numero de colunas.
#
# Tabelas Illumina (8 colunas, sem cabecalho, separadas por tab):
#   bin, chrom, chromStart, chromEnd, name, score, strand, observed
#   - rs: coluna name (5a)
#
# Tabelas Affymetrix (9 colunas, sem cabecalho, separadas por tab):
#   bin, chrom, chromStart, chromEnd, name, score, strand, observed, rsId
#   - rs: coluna rsId (9a)
#   - name contem o identificador interno da sonda Affymetrix
#     (ex.: SNP_A-1886933), nao o rs
#
# Em ambas:
#   - cromossoma: coluna chrom
#   - posicao: coluna chromStart, convertida de zero-based (convencao UCSC)
#     para one-based (convencao do projeto, somando 1)
#   - build: sempre hg19 (descarregados da pasta hg19 do UCSC)
# ---------------------------------------------------------------------------
processar_ucsc <- function(ficheiro, plataforma) {
  
  dados <- read.table(ficheiro,
                      sep = "\t",
                      header = FALSE,
                      quote = "",
                      stringsAsFactors = FALSE,
                      comment.char = "")
  
  if (ncol(dados) == 8) {
    colnames(dados) <- c("bin", "chrom", "chromStart", "chromEnd",
                         "name", "score", "strand", "observed")
    rs_col <- dados$name
    cat("Detetada tabela Illumina (8 colunas)\n")
    
  } else if (ncol(dados) == 9) {
    colnames(dados) <- c("bin", "chrom", "chromStart", "chromEnd",
                         "name", "score", "strand", "observed", "rsId")
    rs_col <- dados$rsId
    cat("Detetada tabela Affymetrix (9 colunas)\n")
    
  } else {
    stop("Estrutura nao reconhecida: ", ncol(dados), " colunas")
  }
  
  # Conversao de zero-based (UCSC) para one-based (convencao do projeto)
  posicao_one_based <- dados$chromStart + 1L
  
  padrao <- criar_padrao(
    rs_bruto       = rs_col,
    chr_bruto      = dados$chrom,
    pos_bruto      = posicao_one_based,
    build_bruto    = "hg19",
    plataforma     = plataforma,
    observed_bruto = dados$observed
  )
  
  gravar_padrao(padrao, plataforma)
  invisible(padrao)
}

# ---------------------------------------------------------------------------
# 4.2. MANIFESTOS ILLUMINA (ficheiros .csv)
# ---------------------------------------------------------------------------
# Estrutura: tem cabecalho mas tambem linhas de metadados no topo, com
# secoes do tipo [Heading], [Assay], [Controls]. Precisamos de localizar
# a linha que comeca com IlmnID para saber onde comecam mesmo os dados.
# Colunas que interessam: Name, Chr, MapInfo, GenomeBuild
#
# Nota: esta funcao ainda nao passa a coluna SNP do manifesto a criar_padrao
# para filtrar indels. A regra de exclusao sera adaptada quando se processarem
# os manifestos Illumina e se confirmar o formato exato da informacao alelica.
# ---------------------------------------------------------------------------
processar_illumina <- function(ficheiro, plataforma) {
  
  # Ler as primeiras 30 linhas para localizar o cabecalho real
  primeiras <- readLines(ficheiro, n = 30)
  linha_cabecalho <- grep("^IlmnID", primeiras)
  
  if (length(linha_cabecalho) == 0) {
    stop("Nao encontrei a linha que comeca por IlmnID nas primeiras 30 linhas")
  }
  
  # Ler o ficheiro a partir do cabecalho verdadeiro
  dados <- read.csv(ficheiro,
                    skip = linha_cabecalho - 1,
                    stringsAsFactors = FALSE,
                    na.strings = c("", "NA", "."),
                    check.names = FALSE)
  
  # Cortar na seccao [Controls] se existir
  if ("IlmnID" %in% colnames(dados)) {
    controls_line <- which(grepl("^\\[Controls\\]", dados$IlmnID))
    if (length(controls_line) > 0) {
      dados <- dados[1:(controls_line[1] - 1), ]
      cat("Encontrada seccao [Controls] na linha", controls_line[1],
          "- dados cortados nessa posicao\n")
    }
  }
  
  cat("Linhas a processar:", nrow(dados), "\n")
  
  # Converter o cromossoma "XY" (regiao pseudoautosomica) para X
  chr_convertido <- as.character(dados$Chr)
  chr_convertido[chr_convertido == "XY"] <- "X"
  
  # Converter a coluna SNP do formato [A/G] para A/G
  snp_limpo <- gsub("[\\[\\]]", "", dados$SNP, perl = TRUE)
  
  # Extrair rs de prefixos compostos como "exm-rs12345" -> "rs12345"
  # Isto recupera rs verdadeiros que estao escondidos atras de prefixos
  # proprietarios dos manifestos Illumina mais antigos (HumanExome, etc.)
  name_limpo <- dados$Name
  padrao_rs_embutido <- "^[a-zA-Z]+-(rs[0-9]+)$"
  com_rs_embutido <- grepl(padrao_rs_embutido, name_limpo)
  if (any(com_rs_embutido)) {
    name_limpo[com_rs_embutido] <- sub(padrao_rs_embutido, "\\1",
                                       name_limpo[com_rs_embutido])
    cat("Identificadores com rs embutido extraidos:",
        sum(com_rs_embutido), "\n")
  }
  
  # Normalizar GenomeBuild: aceitar variantes como "37.1" como "37"
  build_normalizado <- as.character(dados$GenomeBuild)
  build_normalizado <- sub("^37\\..*", "37", build_normalizado)
  
  # Diagnostico antes de chamar criar_padrao
  cat("\nValidacao antes do processamento:\n")
  cat("  Cromossomas distintos apos conversao:\n")
  print(table(chr_convertido, useNA = "ifany"))
  cat("\n  Distribuicao da coluna SNP (10 valores mais frequentes):\n")
  print(head(sort(table(snp_limpo), decreasing = TRUE), 10))
  cat("\n  Distribuicao da GenomeBuild:\n")
  print(table(dados$GenomeBuild, useNA = "ifany"))
  cat("\n")
  
  padrao <- criar_padrao(
    rs_bruto       = name_limpo,
    chr_bruto      = chr_convertido,
    pos_bruto      = dados$MapInfo,
    build_bruto    = build_normalizado,
    plataforma     = plataforma,
    observed_bruto = snp_limpo
  )
  
  # Preservar o nome original (antes da extracao do rs) na coluna nome_original
  # para rastreabilidade
  # Cria correspondencia entre o nome_original do padrao e o Name original
  # Note: esta substituicao funciona porque criar_padrao preserva a ordem,
  # mas com a filtragem alguns indices mudam. Para simplicidade, deixamos
  # nome_original ja com o rs limpo.
  
  gravar_padrao(padrao, plataforma)
  invisible(padrao)
}

# ===========================================================================
# Funcao: processar_illumina_com_aux
# Variante de processar_illumina que integra um ficheiro auxiliar de
# mapeamento Name -> RsID, fornecido pela Illumina para alguns arrays
# (CoreExome v1.0 e v1.1, e potencialmente outros).
#
# A funcao aplica o ficheiro auxiliar APENAS para recuperar rs canonicos
# de marcadores cujo Name original NAO comeca por rs. Nao se aplica
# harmonizacao entre versoes do dbSNP (rs antigo para rs novo), em
# coerencia com a decisao metodologica de deferir essa harmonizacao
# para um passo posterior dedicado.
#
# Regras de aplicacao do ficheiro auxiliar:
#   1. Aplica-se apenas a Names que NAO comecam por rs.
#   2. Aceita-se a recuperacao apenas quando o RsID e um unico rs valido.
#   3. Rejeita-se a recuperacao para RsID descontinuado (".").
#   4. Rejeita-se a recuperacao para RsID com multiplos rs (com virgula).
#   5. Os marcadores rejeitados nao sao perdidos, entram pela via posicional
#      mantendo o Name original.
# ===========================================================================
processar_illumina_com_aux <- function(ficheiro, ficheiro_aux, plataforma) {
  
  # Ler as primeiras 30 linhas para localizar o cabecalho real
  primeiras <- readLines(ficheiro, n = 30)
  linha_cabecalho <- grep("^IlmnID", primeiras)
  
  if (length(linha_cabecalho) == 0) {
    stop("Nao encontrei a linha que comeca por IlmnID nas primeiras 30 linhas")
  }
  
  # Ler o manifesto a partir do cabecalho verdadeiro
  dados <- read.csv(ficheiro,
                    skip = linha_cabecalho - 1,
                    stringsAsFactors = FALSE,
                    na.strings = c("", "NA", "."),
                    check.names = FALSE)
  
  # Cortar na seccao [Controls] se existir
  if ("IlmnID" %in% colnames(dados)) {
    controls_line <- which(grepl("^\\[Controls\\]", dados$IlmnID))
    if (length(controls_line) > 0) {
      dados <- dados[1:(controls_line[1] - 1), ]
      cat("Encontrada seccao [Controls] na linha", controls_line[1],
          "- dados cortados nessa posicao\n")
    }
  }
  
  cat("Linhas a processar:", nrow(dados), "\n")
  
  # Ler o ficheiro auxiliar
  aux <- read.table(ficheiro_aux,
                    sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE,
                    na.strings = c("", "NA"))
  
  cat("Linhas no ficheiro auxiliar:", nrow(aux), "\n")
  
  # Aplicar regras de filtragem ao auxiliar:
  # Aceitar apenas mapeamentos validos para Names nao-rs
  name_nao_rs <- !grepl("^rs", aux$Name)
  rsid_e_rs_unico <- grepl("^rs[0-9]+$", aux$RsID)
  aux_aplicavel <- name_nao_rs & rsid_e_rs_unico
  
  cat("Mapeamentos aplicaveis (Name nao-rs com 1 rs valido):",
      sum(aux_aplicavel), "\n")
  
  # Criar tabela de mapeamento (Name -> rs) apenas com os aplicaveis
  aux_validos <- aux[aux_aplicavel, c("Name", "RsID")]
  
  # Aplicar o mapeamento ao Name do manifesto
  name_original <- dados$Name
  name_atualizado <- name_original
  
  # Fazer o match para substituicao
  posicoes_match <- match(name_original, aux_validos$Name)
  com_mapeamento <- !is.na(posicoes_match)
  name_atualizado[com_mapeamento] <- aux_validos$RsID[posicoes_match[com_mapeamento]]
  
  cat("Names atualizados via auxiliar:", sum(com_mapeamento), "\n")
  
  # Converter o cromossoma "XY" (regiao pseudoautosomica) para X
  chr_convertido <- as.character(dados$Chr)
  chr_convertido[chr_convertido == "XY"] <- "X"
  
  # Converter a coluna SNP do formato [A/G] para A/G
  snp_limpo <- gsub("[\\[\\]]", "", dados$SNP, perl = TRUE)
  
  # Extrair rs embutidos do tipo "exm-rs<numero>" (alem do auxiliar)
  # Isto cobre os casos em que o Name original ja tem rs visivel atras
  # de um prefixo, mas que o auxiliar nao mapeou
  padrao_rs_embutido <- "^[a-zA-Z]+-(rs[0-9]+)$"
  com_rs_embutido <- grepl(padrao_rs_embutido, name_atualizado)
  if (any(com_rs_embutido)) {
    name_atualizado[com_rs_embutido] <- sub(padrao_rs_embutido, "\\1",
                                            name_atualizado[com_rs_embutido])
    cat("Identificadores com rs embutido extraidos (alem do auxiliar):",
        sum(com_rs_embutido), "\n")
  }
  
  # Diagnostico antes de chamar criar_padrao
  cat("\nValidacao antes do processamento:\n")
  cat("  Cromossomas distintos apos conversao:\n")
  print(table(chr_convertido, useNA = "ifany"))
  cat("\n  Distribuicao da coluna SNP (10 valores mais frequentes):\n")
  print(head(sort(table(snp_limpo), decreasing = TRUE), 10))
  cat("\n  Distribuicao da GenomeBuild:\n")
  print(table(dados$GenomeBuild, useNA = "ifany"))
  cat("\n  Names comecados por 'rs' apos atualizacao:",
      sum(grepl("^rs", name_atualizado)), "\n\n")
  
  padrao <- criar_padrao(
    rs_bruto       = name_atualizado,
    chr_bruto      = chr_convertido,
    pos_bruto      = dados$MapInfo,
    build_bruto    = dados$GenomeBuild,
    plataforma     = plataforma,
    observed_bruto = snp_limpo
  )
  
  gravar_padrao(padrao, plataforma)
  invisible(padrao)
}

# ===========================================================================
# Funcao: processar_illumina_anot
# Variante de processar_illumina para Annotation Files da Illumina que
# acompanham arrays sem manifesto CSV padrao acessivel.
#
# Os Annotation Files tem estrutura simples com colunas Name, Chr, MapInfo,
# Alleles, e colunas biologicas adicionais. Nao tem coluna GenomeBuild
# explicita; assume-se hg19 (build 37) por coerencia com a fonte.
#
# Usado para arrays como HumanOmniExpress-24 v1.1 e PsychArray, onde o
# manifesto CSV padrao nao estava disponivel mas o Annotation File estava.
# ===========================================================================
processar_illumina_anot <- function(ficheiro, plataforma) {
  
  # Ler o ficheiro de anotacao
  dados <- read.table(ficheiro,
                      sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE,
                      na.strings = c("", "NA", "."),
                      quote = "", comment.char = "",
                      check.names = FALSE,
                      fill = TRUE)
  
  cat("Linhas a processar:", nrow(dados), "\n")
  cat("Colunas no ficheiro:", paste(colnames(dados), collapse = ", "), "\n\n")
  
  # Verificar se as colunas essenciais estao presentes
  colunas_essenciais <- c("Name", "Chr", "MapInfo", "Alleles")
  colunas_em_falta <- setdiff(colunas_essenciais, colnames(dados))
  if (length(colunas_em_falta) > 0) {
    stop("Colunas essenciais em falta: ",
         paste(colunas_em_falta, collapse = ", "))
  }
  
  # Converter o cromossoma "XY" (regiao pseudoautosomica) para X
  chr_convertido <- as.character(dados$Chr)
  chr_convertido[chr_convertido == "XY"] <- "X"
  
  # Converter a coluna Alleles do formato [A/G] para A/G
  alleles_limpo <- gsub("[\\[\\]]", "", dados$Alleles, perl = TRUE)
  
  # Extrair rs embutidos do tipo "exm-rs<numero>"
  name_limpo <- dados$Name
  padrao_rs_embutido <- "^[a-zA-Z]+-(rs[0-9]+)$"
  com_rs_embutido <- grepl(padrao_rs_embutido, name_limpo)
  if (any(com_rs_embutido)) {
    name_limpo[com_rs_embutido] <- sub(padrao_rs_embutido, "\\1",
                                       name_limpo[com_rs_embutido])
    cat("Identificadores com rs embutido extraidos:",
        sum(com_rs_embutido), "\n")
  }
  
  # Construir coluna build (assume hg19 para todos)
  build_assumido <- rep("37", nrow(dados))
  
  # Diagnostico
  cat("\nValidacao antes do processamento:\n")
  cat("  Cromossomas distintos apos conversao:\n")
  print(table(chr_convertido, useNA = "ifany"))
  cat("\n  Distribuicao da coluna Alleles (10 valores mais frequentes):\n")
  print(head(sort(table(alleles_limpo), decreasing = TRUE), 10))
  cat("\n")
  
  padrao <- criar_padrao(
    rs_bruto       = name_limpo,
    chr_bruto      = chr_convertido,
    pos_bruto      = dados$MapInfo,
    build_bruto    = build_assumido,
    plataforma     = plataforma,
    observed_bruto = alleles_limpo
  )
  
  gravar_padrao(padrao, plataforma)
  invisible(padrao)
}

# ===========================================================================
# Funcao: processar_illumina_anot_com_aux
# Variante de processar_illumina_anot que integra um ficheiro auxiliar
# de mapeamento Name -> RsID, fornecido pela Illumina para alguns arrays
# (PsychArray, e outros que combinem Annotation File com mapeamento aux).
#
# Aplica as mesmas regras metodologicas de processar_illumina_com_aux:
#   1. Aplica-se apenas a Names que NAO comecam por rs.
#   2. Aceita-se a recuperacao apenas quando o RsID e um unico rs valido.
#   3. Rejeita-se a recuperacao para RsID descontinuado (".").
#   4. Rejeita-se a recuperacao para RsID com multiplos rs (com virgula).
#   5. Os marcadores rejeitados nao sao perdidos, entram pela via posicional
#      mantendo o Name original.
# ===========================================================================
processar_illumina_anot_com_aux <- function(ficheiro, ficheiro_aux, plataforma) {
  
  # Ler o ficheiro de anotacao
  dados <- read.table(ficheiro,
                      sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE,
                      na.strings = c("", "NA", "."),
                      quote = "", comment.char = "",
                      check.names = FALSE,
                      fill = TRUE)
  
  cat("Linhas a processar:", nrow(dados), "\n")
  cat("Colunas no ficheiro:", paste(colnames(dados), collapse = ", "), "\n\n")
  
  # Verificar colunas essenciais
  colunas_essenciais <- c("Name", "Chr", "MapInfo", "Alleles")
  colunas_em_falta <- setdiff(colunas_essenciais, colnames(dados))
  if (length(colunas_em_falta) > 0) {
    stop("Colunas essenciais em falta: ",
         paste(colunas_em_falta, collapse = ", "))
  }
  
  # Ler o ficheiro auxiliar
  aux <- read.table(ficheiro_aux,
                    sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE,
                    na.strings = c("", "NA"))
  
  cat("Linhas no ficheiro auxiliar:", nrow(aux), "\n")
  
  # Aplicar regras de filtragem ao auxiliar
  name_nao_rs <- !grepl("^rs", aux$Name)
  rsid_e_rs_unico <- grepl("^rs[0-9]+$", aux$RsID)
  aux_aplicavel <- name_nao_rs & rsid_e_rs_unico
  
  cat("Mapeamentos aplicaveis (Name nao-rs com 1 rs valido):",
      sum(aux_aplicavel), "\n")
  
  aux_validos <- aux[aux_aplicavel, c("Name", "RsID")]
  
  # Aplicar o mapeamento ao Name do annotation
  name_atualizado <- dados$Name
  posicoes_match <- match(name_atualizado, aux_validos$Name)
  com_mapeamento <- !is.na(posicoes_match)
  name_atualizado[com_mapeamento] <- aux_validos$RsID[posicoes_match[com_mapeamento]]
  
  cat("Names atualizados via auxiliar:", sum(com_mapeamento), "\n")
  
  # Converter cromossoma XY para X
  chr_convertido <- as.character(dados$Chr)
  chr_convertido[chr_convertido == "XY"] <- "X"
  
  # Converter Alleles do formato [A/G] para A/G
  alleles_limpo <- gsub("[\\[\\]]", "", dados$Alleles, perl = TRUE)
  
  # Extrair rs embutidos do tipo "exm-rs<numero>"
  padrao_rs_embutido <- "^[a-zA-Z]+-(rs[0-9]+)$"
  com_rs_embutido <- grepl(padrao_rs_embutido, name_atualizado)
  if (any(com_rs_embutido)) {
    name_atualizado[com_rs_embutido] <- sub(padrao_rs_embutido, "\\1",
                                            name_atualizado[com_rs_embutido])
    cat("Identificadores com rs embutido extraidos (alem do auxiliar):",
        sum(com_rs_embutido), "\n")
  }
  
  # Assumir build hg19
  build_assumido <- rep("37", nrow(dados))
  
  # Diagnostico
  cat("\nValidacao antes do processamento:\n")
  cat("  Cromossomas distintos apos conversao:\n")
  print(table(chr_convertido, useNA = "ifany"))
  cat("\n  Distribuicao da coluna Alleles (10 valores mais frequentes):\n")
  print(head(sort(table(alleles_limpo), decreasing = TRUE), 10))
  cat("\n  Names comecados por 'rs' apos atualizacao:",
      sum(grepl("^rs", name_atualizado)), "\n\n")
  
  padrao <- criar_padrao(
    rs_bruto       = name_atualizado,
    chr_bruto      = chr_convertido,
    pos_bruto      = dados$MapInfo,
    build_bruto    = build_assumido,
    plataforma     = plataforma,
    observed_bruto = alleles_limpo
  )
  
  gravar_padrao(padrao, plataforma)
  invisible(padrao)
}

# ===========================================================================
# Funcao: processar_bed
# Processa um ficheiro BED (UCSC/Illumina format) com 4 colunas:
# chr, start (zero-based), end (one-based equivalent), name.
#
# Usada para arrays sem manifesto CSV nem Annotation File acessivel,
# em que so o BED foi possivel obter (caso do Human610-Quad apos
# falha do FTP da Illumina).
#
# LIMITACAO METODOLOGICA IMPORTANTE: o BED nao tem coluna de alelos,
# pelo que NAO e possivel filtrar indels com a regra A/C/G/T usada
# noutros arrays. Esta limitacao e documentada como restricao do
# processamento deste array.
#
# Assume hg19 como build de referencia.
# ===========================================================================
processar_bed <- function(ficheiro, plataforma, ignorar_linhas_track = TRUE) {
  
  # Ler o ficheiro inteiro, ignorando linhas de cabecalho de track
  todas_linhas <- readLines(ficheiro)
  
  # Identificar linhas de track header (comecam por "track" ou "browser")
  if (ignorar_linhas_track) {
    e_track <- grepl("^(track|browser|#)", todas_linhas, ignore.case = FALSE)
    cat("Linhas de cabecalho a ignorar:", sum(e_track), "\n")
    todas_linhas <- todas_linhas[!e_track]
  }
  
  # Dividir cada linha pelos tabs
  campos <- strsplit(todas_linhas, "\t")
  n_campos <- sapply(campos, length)
  
  if (any(n_campos < 4)) {
    cat("AVISO: linhas com menos de 4 campos detectadas, foram ignoradas\n")
    campos <- campos[n_campos >= 4]
  }
  
  # Construir data frame
  dados <- data.frame(
    chr   = sapply(campos, `[`, 1),
    start = as.integer(sapply(campos, `[`, 2)),
    end   = as.integer(sapply(campos, `[`, 3)),
    name  = sapply(campos, `[`, 4),
    stringsAsFactors = FALSE
  )
  
  cat("Linhas a processar:", nrow(dados), "\n")
  
  # Remover prefixo "chr" dos cromossomas (chr1 -> 1, chrX -> X, chrM -> MT)
  chr_limpo <- sub("^chr", "", dados$chr)
  chr_limpo[chr_limpo == "M"] <- "MT"
  chr_limpo[chr_limpo == "XY"] <- "X"
  
  # Converter de zero-based para one-based: posicao = start + 1
  # (equivalente a usar end diretamente para SNPs, mas explicito)
  pos_one_based <- dados$start + 1L
  
  # Extrair rs embutidos do tipo "exm-rs<numero>"
  name_limpo <- dados$name
  padrao_rs_embutido <- "^[a-zA-Z]+-(rs[0-9]+)$"
  com_rs_embutido <- grepl(padrao_rs_embutido, name_limpo)
  if (any(com_rs_embutido)) {
    name_limpo[com_rs_embutido] <- sub(padrao_rs_embutido, "\\1",
                                       name_limpo[com_rs_embutido])
    cat("Identificadores com rs embutido extraidos:",
        sum(com_rs_embutido), "\n")
  }
  
  # Assumir hg19 (sem capacidade de verificar)
  build_assumido <- rep("37", nrow(dados))
  
  # Para arrays sem coluna de alelos, fingir A/G como observed valido
  # para todos os marcadores. Isto desativa o filtro de indels.
  # LIMITACAO: indels (se existirem) NAO serao filtrados.
  observed_fictico <- rep("A/G", nrow(dados))
  
  cat("\nLIMITACAO METODOLOGICA: este array nao tem coluna de alelos\n")
  cat("Indels (se existirem) NAO serao filtrados\n")
  
  # Diagnostico
  cat("\nValidacao antes do processamento:\n")
  cat("  Cromossomas distintos apos conversao:\n")
  print(table(chr_limpo, useNA = "ifany"))
  cat("\n  Primeiros caracteres dos identificadores (top 15):\n")
  print(head(sort(table(substr(name_limpo, 1, 4)), decreasing = TRUE), 15))
  cat("\n")
  
  padrao <- criar_padrao(
    rs_bruto       = name_limpo,
    chr_bruto      = chr_limpo,
    pos_bruto      = pos_one_based,
    build_bruto    = build_assumido,
    plataforma     = plataforma,
    observed_bruto = observed_fictico
  )
  
  gravar_padrao(padrao, plataforma)
  invisible(padrao)
}

# ---------------------------------------------------------------------------
# 4.3. ANOTACOES AXIOM / THERMO FISHER (ficheiros SQLite)
# ---------------------------------------------------------------------------
# Os ficheiros de anotacao dos arrays Axiom da Thermo Fisher vem como bases
# de dados SQLite (extensao .db ou .sqlite), com varias tabelas. A tabela
# de anotacao principal costuma chamar-se Annotations.
#
# Colunas usadas para o pipeline:
#   - dbSNP_RS_ID: identificador rs do dbSNP (NULL para marcadores sem rs)
#   - Chr_id: cromossoma. Convencao Affymetrix: 1-22 nucleares, 24=X, 25=Y,
#     26=MT. Valor 2147483648 (INT_MAX) indica marcadores sem mapeamento.
#   - Start: posicao (convencao one-based, ja em coerencia com o projeto)
#   - Allele_A e Allele_B: alelos. Combinados em pseudo-coluna 'observed'
#     no formato A/B para reutilizar o filtro de criar_padrao. Indels tem
#     '-' num dos alelos e sao automaticamente excluidos.
#
# Conversoes aplicadas:
#   - Cromossomas 24, 25, 26 sao convertidos para X, Y, MT
#   - Build confirmada a partir dos metadados (tabela Information)
#
# Validacoes:
#   - Verifica se a base de dados tem tabelas (deteta ficheiros fantasma)
#   - Compara a build declarada no ficheiro com a indicada como parametro
# ---------------------------------------------------------------------------
processar_axiom <- function(ficheiro_db, plataforma, build = "hg19",
                            tabela = "Annotations") {
  
  con <- dbConnect(RSQLite::SQLite(), ficheiro_db)
  
  # Validar que a base de dados tem tabelas (nao e ficheiro fantasma)
  tabelas <- dbListTables(con)
  if (length(tabelas) == 0) {
    dbDisconnect(con)
    stop("A base de dados nao contem tabelas. Verifica o caminho do ficheiro.")
  }
  
  if (!(tabela %in% tabelas)) {
    dbDisconnect(con)
    stop("Tabela '", tabela, "' nao encontrada. Tabelas disponiveis: ",
         paste(tabelas, collapse = ", "))
  }
  
  # Confirmar a build a partir dos metadados
  if ("Information" %in% tabelas) {
    info <- dbGetQuery(con, "SELECT value FROM Information 
                              WHERE key = 'genome-version-ucsc'")
    if (nrow(info) > 0) {
      cat("Build declarada no ficheiro:", info$value, "\n")
      if (info$value != build) {
        cat("AVISO: build declarada (", info$value,
            ") difere da indicada (", build, ")\n", sep = "")
      }
    }
  }
  
  # Ler colunas relevantes
  dados <- dbGetQuery(con, paste0(
    "SELECT dbSNP_RS_ID, Chr_id, Start, Allele_A, Allele_B FROM ", tabela
  ))
  dbDisconnect(con)
  
  # Forcar conversao para texto antes de qualquer substituicao
  chr_convertido <- as.character(dados$Chr_id)
  
  # Converter cromossomas da convencao Affymetrix para a convencao do projeto
  chr_convertido[chr_convertido == "24"] <- "X"
  chr_convertido[chr_convertido == "25"] <- "Y"
  chr_convertido[chr_convertido == "26"] <- "MT"
  # Valores "2147483648" (INT_MAX) sao marcadores sem mapeamento;
  # serao filtrados pela validacao de cromossoma em criar_padrao
  
  # Construir pseudo-coluna 'observed' no formato A/B esperado por criar_padrao
  observed_construido <- paste(dados$Allele_A, dados$Allele_B, sep = "/")
  
  # Diagnostico antes de chamar criar_padrao
  cat("\nValidacao antes do processamento:\n")
  cat("  Cromossomas distintos apos conversao:\n")
  print(table(chr_convertido, useNA = "ifany"))
  cat("\n  Distribuicao do observed construido (10 valores mais frequentes):\n")
  print(head(sort(table(observed_construido), decreasing = TRUE), 10))
  cat("\n")
  
  padrao <- criar_padrao(
    rs_bruto       = dados$dbSNP_RS_ID,
    chr_bruto      = chr_convertido,
    pos_bruto      = dados$Start,
    build_bruto    = build,
    plataforma     = plataforma,
    observed_bruto = observed_construido
  )
  
  gravar_padrao(padrao, plataforma)
  invisible(padrao)
}

# ===========================================================================
# 5. JUNCAO FINAL DOS FICHEIROS PADRAO
# ===========================================================================
# Quando tiver todos os ficheiros padrao_*.tsv prontos, este bloco junta
# tudo num unico ficheiro mestre para o cruzamento posterior.
# ===========================================================================
juntar_padroes <- function(pasta = ".", saida = "padrao_TODOS.tsv") {
  
  ficheiros <- list.files(pasta,
                          pattern = "^padrao_.*\\.tsv$",
                          full.names = TRUE)
  
  ficheiros <- ficheiros[!grepl(saida, ficheiros)]
  
  cat("Ficheiros a juntar:", length(ficheiros), "\n")
  
  mestre <- do.call(rbind, lapply(ficheiros, function(f) {
    read.table(f, sep = "\t", header = TRUE,
               stringsAsFactors = FALSE, na.strings = "NA",
               quote = "", comment.char = "")
  }))
  
  write.table(mestre, file = saida,
              sep = "\t", row.names = FALSE,
              quote = FALSE, na = "NA",
              fileEncoding = "UTF-8")
  
  cat("Ficheiro mestre gravado em:", saida, "\n")
  cat("Total de linhas:", nrow(mestre), "\n")
  cat("Plataformas incluidas:\n")
  print(table(mestre$plataforma))
  
  invisible(mestre)
}