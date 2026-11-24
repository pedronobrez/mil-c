---
title: Histórico de versões
---

## 1.0.0 — 24 de novembro de 2026

A primeira versão. Ela cobre o caminho que veio substituir — gatear no FlowJo, DownSample,
concatenar, rodar os plugins de UMAP, FlowSOM e ClusterExplorer — e vai bem além dele. O que a
1.0 contém, por área:

### Arquivos e amostras

- FCS 2.0, 3.0 e 3.1 lidos e gravados só com numpy; workspaces do FlowJo 10 lidos — gates,
  compensação e as contagens do FlowJo para comparar — inclusive quadrantes, combinações
  booleanas, elipses e quadrantes curvos ([[opening-data]]).
- A tabela da aba Samples é sua para digitar: colunas livres, valores digitados em várias linhas
  de uma vez, colunas lidas dos nomes dos arquivos com um padrão, um grupo feito de uma coluna ou
  de várias ([[samples-tab]]).
- Uma ordem dos grupos para todo gráfico, legenda, painel, tabela e arquivo do Prism, guardada com
  o projeto ([[sample-order]]).
- Tudo o que se faz com uma pasta guardado ao lado dos arquivos em `MIL-C.json`
  ([[projects-and-saving]]).
- Tubos readquiridos apontados antes de contarem em dobro.

### Concatenação, normalização, qualidade

- Concatenação balanceada com o plano na mesa antes de qualquer sorteio: o mesmo número de cada
  amostra, o mesmo total de cada grupo, ou um teto; um mínimo; um sorteio com semente que não
  depende das outras amostras ([[concatenation]]).
- Normalização de lotes à maneira do CytoNorm, em cada marcador inteiro ou dentro de clusters do
  FlowSOM, para a média dos lotes ou para um lote, opcionalmente aplicada às próprias amostras;
  uma vista de antes e depois sobre eventos separados do ajuste ([[normalisation]]).
- Controle de qualidade de cada aquisição ao longo do tempo — vazão, sinal, derivas que duram e
  desvios curtos — como uma população Clean no topo de cada árvore, com os limites no diálogo
  ([[quality-control]]).

### Gates

- Gates de retângulo, polígono, laço, elipse, faixa e quadrante, quadrantes cujos braços se movem
  um a um, combinações booleanas, quadrantes curvos desenhados e salvos; um gate para todas as
  amostras ou para uma só ([[drawing-gates]]).
- Gates automáticos — um vale, um contorno de densidade, uma elipse robusta, uma faixa de
  singletos — reajustados sozinhos quando a população em que estão muda ([[auto-gates]]).
- Escalas biexponencial, Logicle, arcsinh, logarítmica e linear por parâmetro
  ([[scales-and-transforms]]).
- Uma matriz de compensação por amostra ou comum, editada numa tabela ([[compensation]]); uma
  matriz a partir de controles de uma cor pelo método do AutoSpill ou pelo clássico
  ([[autocompensation]]); unmixing espectral por OLS ou WLS ([[spectral-unmixing]]).
- Sobreposições, violinos e caixas, uma matriz com todo par de parâmetros, backgating
  ([[gates-tab]]).
- Modelos de uma população: ciclo celular ([[cell-cycle]]), proliferação ([[proliferation]]) e
  cinética ([[kinetics]]).
- Gates exportados como workspace do FlowJo ou como Gating-ML 2.0
  ([[workspaces-and-gating-ml]]).

### Estatística

- Uma tabela de estatísticas de toda população de toda amostra: contagens, porcentagens,
  medianas, médias, CVs e medidas robustas pelas definições do FlowJo ([[statistics-tab]]).
- Testes entre grupos, desenhos pareados e de medidas repetidas, ANOVA de Welch, com correção para
  muitos testes.
- Tabelas enviadas ao GraphPad Prism como arquivos `.pzfx` ([[tables-and-prism]]).
- Abundância diferencial dos clusters à maneira do DA-edgeR do diffcyt ([[abundance]]) e
  expressão diferencial dentro dos clusters à maneira do DS-limma
  ([[differential-expression]]).

### Explore e Clusters

- UMAP, tSNE, FlowSOM, agrupamento no estilo do PhenoGraph com Leiden ou Louvain e tamanho mínimo
  de cluster, e PCA, rodados em segundo plano ([[embeddings]], [[clustering-methods]]).
- Mapas coloridos e divididos por qualquer coluna, ou por duas ([[explore-tab]]).
- O mapa, o heatmap, a árvore e as frequências dos clusters, ligados; clusters nomeados,
  anotados, juntados e transformados em populações; testes por cluster ([[clusters-tab]],
  [[cluster-statistics]]).

### Figuras e relatórios

- Todo gráfico como figura em milímetros, com fonte, cores e rótulos, em SVG, PDF, PNG ou TIFF;
  figuras da aba Gates levam os seus gates ([[figures]]).
- Pranchas de vários painéis, salvas com o projeto, e um relatório em PDF de como a análise foi
  feita ([[the-plate]]).

### Plataforma

- Instaladores para macOS (Apple Silicon), Windows (`.msi` e um `.zip` portátil) e Linux (`.deb`
  e `.tar.gz`), montados e testados por integração contínua nos três sistemas, com checksums
  ([[installation]]).
- Uma linha de comando que faz o que a janela faz, e um autoteste ([[command-line]]).
- Aparência clara e escura; um layout que cabe na tela de um notebook e nas molduras mais largas
  do Windows.
- Este manual, em inglês e em português, aberto na página de onde você está, e exportável como um
  único PDF ([[welcome]]).

### Conferido com

O FlowJo em aquisições reais, e os pacotes de R e Python que cada método segue — edgeR,
PhenoGraph, PeacoQC, flowAI, CytoNorm, AutoSpill, CytoML, os valores exportados pelo SpectroFlo e
o AutoSpectral. Os resultados, com onde cada um foi medido, estão em [[measured-facts]].
