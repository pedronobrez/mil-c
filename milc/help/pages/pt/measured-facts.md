---
title: Medido contra outros programas
---

A maior parte dos métodos do MIL-C foi escrita a partir dos artigos, não copiada dos programas que os tornaram conhecidos. Cada um foi depois posto ao lado do programa em que se inspira, com a mesma entrada, e o que está escrito aqui é o que foi medido. Onde um número não foi medido, ele não aparece; onde algo não foi conferido, esta página diz. O raciocínio por trás disso está em [[design-principles]].

As comparações com pacotes de R e Python foram rodadas em 23 de novembro de 2026, com R 4.6.1, Bioconductor 3.23 e Python 3.14. Os scripts estão em `tests/reference/` no código-fonte; eles precisam de R e não fazem parte da suíte de testes. "Um experimento real", abaixo, quer dizer aquisições que não são distribuídas com o MIL-C; os conjuntos de dados públicos são nomeados.

## Contra o FlowJo

Num experimento real de 12 amostras com gates e cinco marcadores, e dois arquivos concatenados de cerca de 2.000 eventos cada. A suíte de testes confere estes números com uma margem.

| O quê | MIL-C contra o FlowJo |
|---|---|
| Contagens dos gates, 72 gates | Dentro de 0,77% em todos os gates; dentro de 3 eventos na população final. O FlowJo avalia polígonos numa grade de 256 posições por uma regra que não foi reproduzida; o MIL-C usa a geometria exata |
| Transformação biexponencial | Medianas por cluster dentro de 0,69 de um canal (de 1.023) |
| Clusters do FlowSOM (índice de Rand ajustado) | 0,86 e 0,81 contra os do FlowJo; duas execuções do MIL-C com sementes diferentes concordam em 0,81 e 0,79. A execução do FlowJo fica tão perto da do MIL-C quanto duas do MIL-C entre si |
| UMAP: fração dos 15 vizinhos mais próximos de cada evento mantida | 47,0% e 47,9%; FlowJo 47,2% e 47,7% |
| tSNE: a mesma medida | 58,8%; FlowJo 58,0% |

O UMAP é o `umap-learn`, não o `uwot` que o plugin do FlowJo chama: as imagens diferem, as vizinhanças que elas mantêm não. Ver [[embeddings]] e [[clustering-methods]].

**Não postos ao lado do FlowJo:** as [[statistics-tab|estatísticas]] (seguem as definições do FlowJo); os quadrantes "curly" (o FlowJo não grava a curva); proliferação e ciclo celular (nem ao lado do ModFit). Elipses lidas de um workspace foram conferidas uma vez, num arquivo publicado de 100 eventos com eixos lineares.

## Abundância diferencial contra o edgeR e o diffcyt

A [[abundance|abundância diferencial]] ao lado do edgeR 4.10.7 e do diffcyt 1.32.1, em 600 tabelas de contagens simuladas (7.600 clusters: 5 contra 5 e 3 contra 3 amostras, 8 ou 15 clusters, metade com três clusters dobrados) e em quatro tabelas de um experimento real de dois grupos (46 clusters).

| edgeR chamado | Razão (log2) | Valor de p | Descobertas a FDR 0,05 |
|---|---|---|---|
| Como o MIL-C imita (dispersão comum, teste QL anterior ao edgeR 4) | Dentro de 0,0004 | Dentro de 0,06% em 95% dos clusters | 616 em comum, 0 diferentes |
| Como o MIL-C imita (dispersão por cluster, 10 graus a priori, LRT) | Dentro de 0,0003 | Dentro de 2% em 95% | 724 em comum, 1 diferente |
| Padrões do edgeR 4 (QL) | Mediana 0,0005 | Um fator de 2 em 5% dos clusters | Cerca de uma descoberta em cinco difere |
| `diffcyt::testDA_edgeR` (LRT) | Mediana 0,00008 | Um fator de 2 em 5% | Cerca de uma descoberta em dez difere |

A dispersão comum fica dentro de 0,4% da do edgeR. A diferença para os padrões é de método: eles estimam os graus de liberdade a priori das dispersões, o MIL-C os fixa em 10. A comparação achou dois defeitos no teste F do MIL-C (graus de liberdade não limitados aos somados; a variância comum tomada por média geométrica), e os dois foram corrigidos. Calibração nas mesmas tabelas: o edgeR tem o mesmo leve excesso de valores de p pequenos sem diferença que o MIL-C, e nos seus padrões um pouco mais.

## Expressão diferencial contra o limma

A [[differential-expression|expressão diferencial]]: coeficientes, t moderado e p coincidem com os do limma 3.68.5 em cerca de sete algarismos significativos numa matriz 8 × 8 em seis desenhos (com e sem covariável, pesos e valores faltando). Com medianas faltando ela segue `eBayes(legacy = TRUE)`; o estimador mais novo do limma dá graus de liberdade a priori um pouco diferentes (3,52 contra 3,10 na matriz de teste). O diffcyt em si não foi rodado.

## PhenoGraph contra o pacote dos autores

O [[clustering-methods|PhenoGraph]] ao lado do PhenoGraph 1.5.7 (com o seu Louvain e com Leiden), k = 30, resolução 1, três execuções de cada:

- Em misturas de seis gaussianas em 8 dimensões (populações de 40% a 1%): o Leiden do MIL-C e o do pacote deram as mesmas partições nas misturas bem separadas (ARI 1,00 e 0,997); nas sobrepostas, 0,75 a 0,78, onde o pacote concorda consigo mesmo em 0,75 a 0,77.
- Em 20.000 eventos reais com seis marcadores (medido quando o MIL-C ainda usava Louvain): o mesmo número de clusters, com diferença de um ou dois, e o MIL-C concorda com o pacote (ARI de cerca de 0,7) tanto quanto cada um concorda consigo mesmo de uma execução para outra.
- A modularidade do Leiden ao lado da do leidenalg no mesmo grafo: 0,8065 contra 0,8036.
- O MIL-C levou de 1 a 7 segundos por execução, o pacote de 27 a 183 (não é uma medida de desempenho cuidadosa).

## Controle de qualidade contra o PeacoQC e o flowAI

O [[quality-control|controle de qualidade]] ao lado do PeacoQC 1.22.0 e do flowAI 1.42.0:

- **Defeitos conhecidos** (sintéticos e escritos em arquivos reais): um degrau ou uma queda brusca são removidos como o PeacoQC os remove; um entupimento é removido inteiro, onde o PeacoQC, que não examina a vazão, não remove nada; de uma deriva lenta, 76% a 97% são removidos contra 74% a 98% do PeacoQC. Do resto do arquivo, fora do defeito, o MIL-C removeu de 0 a 0,7%, o PeacoQC de 0 a 2%. O flowAI, nos valores brutos e nos seus padrões, achou poucos dos defeitos escritos.
- **32 arquivos reais sem defeito conhecido:** o MIL-C remove 5,5% dos eventos, o PeacoQC 5,1%, o flowAI 12,0% — não os mesmos eventos: do que o PeacoQC remove, o MIL-C remove 56%. Três arquivos dos quais o PeacoQC remove de 12% a 31% ficam intactos, sem que a mediana de nenhum canal se mova ali mais que 0,05 do seu espalhamento. Sem a verdade, isso não diz quem tem razão.
- **Cópias de arquivos reais com um defeito de tamanho conhecido:** um desvio de 0,3 do espalhamento durante 5% do arquivo é achado em 24 de 32; o PeacoQC o achou em 8 de 8 arquivos menores, mas removeu 6,3% dos mesmos arquivos sem defeito, onde o MIL-C remove 0,04%. Os arquivos reais perdem o mesmo com o limite de deriva em qualquer ponto de 0,20 a 0,33 e o de desvio curto de 0,20 a 0,30: o padrão 0,25 fica no meio das duas faixas.

## Normalização de lotes contra o CytoNorm

A [[normalisation|normalização]] ao lado do CytoNorm 2.0.13, em lotes simulados com efeito conhecido e num experimento real de três lotes. Sem clusters, é a curva do CytoNorm dentro da faixa coberta pelas referências (diferença média de cerca de 1e-5 do eixo; um milésimo da correção feita). Com clusters, os dois programas ficam tão longe um do outro quanto cada um de si mesmo com outra semente. Além dos quantis extremos eles diferem de propósito: o CytoNorm continua a inclinação da spline, o MIL-C continua com um deslocamento.

## Compensação por controles contra o AutoSpill

A [[autocompensation|autocompensação]] ao lado do pacote R AutoSpill (0.2.0) e do `flowStats::spillover`: nos mesmos eventos, a matriz do AutoSpill a 2e-5 em controles simulados e a 3e-4 nos controles de uma cor reais que vêm com o flowCore, no mesmo número de iterações. Em controles simulados com matriz conhecida, o MIL-C fica dentro de 0,04 ponto percentual (o método clássico, dentro de 0,3). Nada foi medido num aparelho moderno de vários lasers, em beads ou num espectral, e nenhuma matriz foi posta ao lado da do FlowJo ou do FACSDiva.

## Unmixing espectral contra o SpectroFlo e o AutoSpectral

O [[spectral-unmixing|unmixing]], num conjunto público de beads do Aurora (42 controles, 64 detectores) com todo arquivo bruto e também desmisturado pelo SpectroFlo 3.3: com os espectros do SpectroFlo, o MIL-C reproduz os valores desmisturados dele a 1e-7; os espectros que o MIL-C tira dos mesmos controles brutos têm cosseno mediano de 0,99994 com os do SpectroFlo (o menor, 0,966); as razões das medianas desmisturadas das beads ficam entre 0,991 e 1,023. OLS, WLS por evento, a matriz de similaridade e o número de condição coincidem com as funções R do AutoSpectral a 1e-9 ou melhor. O WLS e o índice de complexidade próprios da Cytek não são publicados e não foram comparados.

## Ciclo celular, proliferação e cinética

- [[cell-cycle|Ciclo celular]], em 223 arquivos públicos de BrdU/PI: ajustado só no PI, os agregados achados nos eventos sem gate batem com o que um gate de singletos PI-A/PI-H retira (24,8% contra 25,6%, r = 0,90). Nos 39 controles, o Dean-Jett-Fox dá 6,6 pontos de S a menos que as células BrdU+ (DP 5,8) — o BrdU também conta células que saíram de S durante o pulso — e o Watson 3,9 pontos a mais (DP 13,5).
- A [[proliferation|proliferação]] coincide com o flowFit do R nos próprios dados dele de CFSE e CellTrace Violet, dentro de 0,5 ponto por geração; na amostra de Cell Proliferation Dye os dois diferem em até 6,5 pontos.
- A [[kinetics|cinética]] coincide com as mesmas definições escritas de novo em R a 1e-13. Isso confere a aritmética, não um fluxo de cálcio real.

## Gates, escalas e arquivos

- A Logicle coincide com a do FlowUtils em 15 casas decimais (ver [[scales-and-transforms]]).
- Quadrantes "curly", num workspace escrito pelo MIL-C com 40.000 eventos simulados, põem os mesmos eventos em cada canto que o CytoML 2.24.0, a menos de até 65 de um canto, todos vindos da tabela mais grossa da biexponencial do cytolib; com aquela tabela no lugar, no máximo 1.
- Um workspace escrito pelo MIL-C foi lido de volta pelo MIL-C e pelo FlowKit e deu os mesmos eventos (32 amostras, 192 populações), e foi aberto à mão no FlowJo antes da 1.0. Os arquivos Gating-ML validam contra o esquema do padrão (ver [[workspaces-and-gating-ml]]).
- [[auto-gates|Gates automáticos]]: nos dados de teste, os gates de vale coincidem com os desenhados à mão com Jaccard de 0,95 a 0,97; os gates de dispersão saem mais folgados.
- Os arquivos do Prism seguem o pacote R `pzfx` e fazem a ida e volta pelo leitor do MIL-C; foram abertos à mão no Prism antes da 1.0 (ver [[tables-and-prism]]).
- A ANOVA de Welch e a de medidas repetidas reproduzem exemplos de livro e do R e cálculos independentes em dados aleatórios (ver [[cluster-statistics]]).
