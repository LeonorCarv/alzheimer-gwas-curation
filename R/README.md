# Scripts em R

Processamento das plataformas de genotipagem e análise de sobreposição de SNPs
entre plataformas e entre GWAS. Os scripts correm em cadeia, cada um a produzir
ficheiros que o seguinte vai ler.

## Ordem de execução

O `pipeline.R` corre tudo de uma vez, pela ordem certa, e é o ponto de entrada
recomendado. Por baixo, a sequência é esta:

1. `funcoes_padronizacao.R` reúne as funções que normalizam os manifestos das
   várias fontes (UCSC, Illumina, Axiom da Thermo Fisher) e os juntam num
   ficheiro mestre único.
2. `analise_cruzamento_platf.R` calcula a sobreposição de SNPs entre plataformas.
3. `preparar_gwas.R` estrutura a tabela de GWAS a partir do CSV de metadados.
4. `analise_gwas.R` calcula, para cada GWAS, os SNPs por união e por interseção
   das plataformas usadas.
5. `analise_gwas_cruzamento.R` constrói as matrizes de sobreposição de SNPs entre
   GWAS, dois a dois.
6. `generate_supplementary_excel.R` reúne tudo num único ficheiro Excel com as
   tabelas suplementares.

## Dados de entrada

Os manifestos das plataformas (UCSC, Illumina, Axiom) e a tabela de metadados dos
GWAS. Os manifestos não estão no repositório por serem ficheiros grandes e de
fontes externas. O cabeçalho do `pipeline.R` lista quais são e onde os obter.

## Como correr

A partir da raiz do projeto, em R:

    source("R/pipeline.R")

Os caminhos dos ficheiros são relativos, por isso convém correr sempre a partir
da mesma pasta.

## Dependências

R com os pacotes *zoo*, *DBI*, *RSQLite* e *openxlsx*.
