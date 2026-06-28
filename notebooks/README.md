# Notebook graphs

Visualizações do projeto de curadoria de GWAS de doença de Alzheimer. O notebook
parte das tabelas de metadados dos estudos e gera os gráficos e a tabela usados
na análise.

## Conteúdo

- `graphs.ipynb` produz, por esta ordem, os gráficos de ligação entre bases de
  dados e GWAS (casos e controlos), os gráficos de ligação entre métodos de
  diagnóstico e GWAS (casos e controlos), o heatmap da sobreposição de indivíduos entre
  estudos, a tabela de bases de dados por metodologia e o gráfico de barras
  empilhado com a cobertura de SNPs por GWAS.

## Dados de entrada

O notebook lê dois ficheiros, ambos separados por ponto-e-vírgula e em
codificação cp1252:

- `metadata_gwas.csv`, a tabela de metadados dos GWAS
- `metadata_snps.csv`, a tabela de SNPs por estudo

Os dois são lidos da pasta de trabalho, sem caminho indicado. Como no
repositório estão em `data/raw/`, para correr a partir desta pasta tem de se 
copiar para junto do notebook ou ajustar o caminho nas células de leitura.

## Ficheiros gerados

Ao correr, o notebook cria as pastas de saída automaticamente:

- `created_figuras/` com as figuras em PNG a 300 dpi, organizadas em subpastas
  por tipo de gráfico (`sobrep_BD_GWAS`, `sobrep_Diag_GWAS`, `snpsVSgwas`)
- `created_tables/` com a tabela de estudos por plataforma, em CSV e HTML

## Como correr

Abrir o notebook no Jupyter e executar as células pela ordem em que aparecem
(menu Run, opção Run All).

## Dependências

Python 3 com *pandas*, *NumPy* e *matplotlib*. O resto é biblioteca padrão.
