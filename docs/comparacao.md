# Quatro métodos do MIL-C ao lado dos pacotes originais

Medido em 23 de novembro de 2026. Quatro métodos do MIL-C foram escritos a partir dos artigos e
só tinham sido conferidos por simulação. Aqui cada um foi posto ao lado do programa em que se
inspira, com a mesma entrada, e o que está escrito é o que foi medido. Onde um número não foi
medido, ele não aparece.

Os scripts que produzem tudo isto estão em `tests/reference/` (o `README.md` de lá diz como
instalar os pacotes e rodar). Eles não fazem parte da suíte de testes.

## Resumo

| Método do MIL-C | Comparado com | Resultado em uma frase |
|---|---|---|
| `abundance.differential` | edgeR 4.10.7, diffcyt 1.32.1 | Chamado do mesmo jeito, o edgeR dá os mesmos números (p dentro de 0,06% no teste QL e 2% no LRT em 95% dos clusters). Contra os *padrões* do edgeR 4 e do diffcyt a diferença é maior, por desenho. A comparação achou dois defeitos no teste QL, que foram corrigidos. |
| `graphcluster.phenograph` | PhenoGraph 1.5.7 (PyPI) | Mesmo número de clusters, com diferença de um ou dois, e partições tão parecidas com as do pacote quanto duas corridas do pacote entre si (ARI 0,71 nos eventos reais; 0,73 pacote contra pacote). |
| `qc.check` | PeacoQC 1.22.0, flowAI 1.42.0 | Defeitos nítidos (degrau, entupimento) são achados igual ou melhor. Desvios lentos e moderados o PeacoQC achava e o nosso não (1,4% dos eventos reais removidos contra 5,1%): `qc.check` ganhou uma passada lenta e, medido de novo, acha-os no nível do PeacoQC (5,25% contra 5,13% dos eventos reais; 76% contra 74% de uma deriva injetada), sem as remoções dele em arquivos cujas medianas não se movem. |
| `normalize.fit` / `apply` | CytoNorm 2.0.13 (GitHub) | Sem clusters é a mesma curva (diferença média de 1e-5 do eixo). Com clusters, os dois programas diferem entre si tanto quanto duas sementes de um mesmo. Fora da faixa das referências o nosso é deliberadamente mais contido. |
| `autocomp.autospill` / `classic` | AutoSpill 0.2.0 (GitHub), flowStats 4.24.0 | Nos mesmos eventos o nosso AutoSpill dá a matriz do pacote a 2e-5 (simulado) e 3e-4 (controles reais do flowCore), com o mesmo número de iterações. Em simulação acha a matriz verdadeira a 0,04 ponto percentual; o clássico, a 0,3. Os exemplos do próprio AutoSpill (FlowRepository) não puderam ser baixados. |
| `cellcycle.fit` | o BrdU das mesmas células (223 arquivos públicos de BrdU/PI) | O FlowJo e o ModFit não estão aqui. Sem gate de singletos, os agregados achados (24,8%) batem com o que o gate tira (25,6%; r = 0,90) e as fases ficam as dos singletos. Nos 39 controles o Dean-Jett-Fox dá 6,6 pontos de S a menos que o BrdU (DP 5,8); o Watson, 3,9 a mais com DP 13,5. |
| `kinetics.curve` / `measure` | as mesmas definições refeitas em R (flowCore) | Iguais até 1e-13 em 23 arquivos reais com Time. Confere a aritmética; nenhum fluxo de cálcio público foi alcançado. |

Versões: R 4.6.1, Bioconductor 3.23, limma 3.68.5, flowCore 2.24.0, FlowSOM 2.20.0, leidenalg
0.12.0, Python 3.14.6, em macOS arm64. Todos os pacotes instalaram; nenhuma comparação foi
pulada.

Dados reais: o experimento de 30 amostras (dias 3, 7 e 14; Case e Ctrl; 5 animais por grupo), com
os gates do workspace do FlowJo. As duas aquisições repetidas (com crase no nome) ficaram fora,
menos no controle de qualidade, que olha arquivo por arquivo e usou as 32.

---

## 1. Abundância diferencial (`milc/abundance.py`) e o edgeR

### O que foi comparado

A mesma tabela de contagens (clusters × amostras) foi dada ao nosso `differential` e ao edgeR, e
comparados: log2 da razão, dispersão, valor de p e FDR (Benjamini-Hochberg), cluster a cluster. O
edgeR foi chamado de cinco maneiras:

| Nome | Chamadas | Ao lado do nosso |
|---|---|---|
| `ql_matched` | `estimateGLMCommonDisp`; `glmQLFit(dispersion = comum, legacy = TRUE, abundance.trend = FALSE)`; `glmQLFTest` | QL |
| `ql_default` | `estimateDisp`; `glmQLFit`; `glmQLFTest`, tudo no padrão do edgeR 4 | QL |
| `lrt_matched` | `estimateGLMCommonDisp`; `estimateGLMTagwiseDisp(prior.df = 10, trend = FALSE)`; `glmFit`; `glmLRT` | LRT |
| `lrt_disp` | `estimateDisp(trend.method = "none")`; `glmFit`; `glmLRT` — as chamadas que o diffcyt faz | LRT |
| `diffcyt` | `diffcyt::testDA_edgeR`, sem filtro de clusters raros | LRT |

As duas `matched` são o que o nosso código diz imitar. As outras três são o que um usuário do
edgeR ou do diffcyt obtém hoje sem mexer em nada. Em todas o tamanho da biblioteca é o total da
amostra, sem fatores de normalização (o padrão do diffcyt).

### Dados

- **Simuladas**: 600 tabelas da simulação de `tests/test_abundance.py` (clusters de 0,1% a 30%,
  dispersão própria de 0,02 a 0,3, totais variando cinco vezes). São 100 sem diferença e 100 com
  três clusters dobrados, em cada um de três cenários: 5 contra 5 com 15 clusters, 3 contra 3
  com 15, 5 contra 5 com 8. Ao todo 7.600 clusters.
- **Reais**: os eventos CD45+ das 30 amostras (cerca de 80 mil), agrupados em 12 metaclusters
  pelo FlowSOM do MIL-C em seis marcadores, contados por amostra. Quatro tabelas: Case contra
  Ctrl em cada dia (5 contra 5) e nos três dias juntos com o dia como covariável (15 contra 15).
  São 46 clusters testados.

### Concordância medida

Tabelas simuladas, nosso contra cada variante do edgeR (diferenças absolutas):

| Variante do edgeR | log2FC: mediana | log2FC: máx. | Dispersão, razão em log2: mediana | p95 | log10 p: mediana | log10 p: p95 | log10 p: máx. | Clusters com p diferente em mais de 0,1 dex | FDR<0,05 nos dois | Só o nosso | Só o edgeR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `ql_matched` | 3,3e-06 | 0,0004 | 0,0042 | 0,0056 | 2,2e-05 | 0,00025 | 0,60 | 2 | 616 | 0 | 0 |
| `lrt_matched` | 3,6e-06 | 0,0003 | 0,0047 | 0,0098 | 0,0005 | 0,0087 | 0,077 | 0 | 724 | 0 | 1 |
| `ql_default` | 0,0005 | 0,16 | 0,041 | 0,22 | 0,027 | 0,37 | 4,4 | 1.650 | 592 | 24 | 120 |
| `lrt_disp` e `diffcyt` | 8e-05 | 0,041 | 0,12 | 0,57 | 0,015 | 0,31 | 21 | 1.118 | 696 | 28 | 46 |

Uma diferença de 0,00025 em log10 p é 0,06% no valor de p; 0,0087 é 2%; 0,3 é um fator de 2.
`lrt_disp` e `diffcyt` deram valores de p idênticos entre si, como esperado.

A dispersão comum do nosso fica 0,3% acima da do `estimateGLMCommonDisp` (mediana; no máximo
0,45%): o nosso procura o máximo numa grade com uma parábola, o edgeR com `optimize`.

Tabelas reais:

| Variante do edgeR | log2FC: mediana | log2FC: máx. | log10 p: mediana | log10 p: p95 | log10 p: máx. | FDR<0,05 nos dois | Só o nosso | Só o edgeR |
|---|---|---|---|---|---|---|---|---|
| `ql_matched` | 8e-06 | 0,003 | 6e-05 | 0,08 | 0,57 | 8 | 0 | 0 |
| `lrt_matched` | 3e-05 | 1,16 | 0,0005 | 0,055 | 0,56 | 11 | 0 | 0 |
| `ql_default` | 0,002 | 0,57 | 0,10 | 1,07 | 5,1 | 8 | 0 | 0 |
| `lrt_disp` e `diffcyt` | 0,001 | 2,7 | 0,07 | 1,5 | 11,7 | 8 | 3 | 0 |

As mesmas 600 tabelas simuladas, cada método por si (o quanto cada um acerta, não o quanto se
parecem):

| Cenário | Método | p<0,05 sob a nula (%) | Poder a FDR<0,05 (%) | FDR obtida (%) |
|---|---|---|---|---|
| 5 v 5, 15 clusters | nosso QL | 4,5 | 69 | 6,0 |
| | edgeR `ql_matched` | 4,5 | 69 | 6,0 |
| | edgeR `ql_default` | 5,3 | 72 | 7,4 |
| | nosso LRT | 5,1 | 73 | 7,6 |
| | edgeR `lrt_matched` | 5,1 | 73 | 7,6 |
| | diffcyt | 5,0 | 74 | 7,7 |
| 3 v 3, 15 clusters | nosso QL | 6,5 | 34 | 4,3 |
| | edgeR `ql_matched` | 6,5 | 34 | 4,3 |
| | edgeR `ql_default` | 6,9 | 42 | 9,0 |
| | nosso LRT | 6,3 | 49 | 9,0 |
| | edgeR `lrt_matched` | 6,3 | 49 | 9,0 |
| | diffcyt | 6,5 | 49 | 9,8 |
| 5 v 5, 8 clusters | nosso QL | 5,9 | 67 | 3,8 |
| | edgeR `ql_matched` | 5,9 | 67 | 3,8 |
| | edgeR `ql_default` | 6,5 | 70 | 5,6 |
| | nosso LRT | 6,4 | 73 | 3,4 |
| | edgeR `lrt_matched` | 6,4 | 73 | 3,4 |
| | diffcyt | 6,8 | 73 | 5,0 |

Ou seja: o excesso de valores de p pequenos que a nossa simulação já mostrava não é um defeito
da nossa implementação. O edgeR tem o mesmo excesso nas mesmas tabelas, e nos seus padrões um
pouco mais.

### Onde diferem e por quê

**Dois defeitos nossos, corrigidos** (em `milc/abundance.py`, com teste em
`tests/test_abundance.py::test_the_numbers_are_edgers`, que guarda duas tabelas e os números que
o edgeR deu para elas):

1. *Graus de liberdade do teste F.* Quando as quase-dispersões dos clusters não diferem mais do
   que o acaso explicaria, a priori do `squeezeVar` vale "infinitos" graus de liberdade, e o
   nosso código testava contra um qui-quadrado. O limma e o edgeR limitam os graus de liberdade
   totais aos graus somados de todos os clusters (clusters × resíduo), porque a variância comum
   foi estimada com isso e não com infinito. Aconteceu em 8% das tabelas simuladas. Na tabela
   real do dia 7 o nosso p do cluster mais alterado era 2,7e-31 e o do edgeR 2,5e-20.
2. *A variância comum nesse mesmo caso.* O nosso usava uma média geométrica corrigida; o limma
   (`fitFDist`) usa a média aritmética. A estatística F diferia em até 20%.

Antes da correção, nas mesmas 600 tabelas, 53 clusters tinham p diferente do `ql_matched` em
mais de 0,1 dex e 4 descobertas a FDR<0,05 não coincidiam; depois, 2 clusters e nenhuma
descoberta. Nas tabelas reais o máximo caiu de 11 dex para 0,57. O efeito prático é pequeno: a
fração de p<0,05, p<0,01 e p<0,001 sob a nula ficou igual antes e depois nessas tabelas.

**Diferenças que ficam, por desenho:**

- *`poisson.bound`.* Onde a variância ajustada fica abaixo da de um Poisson, o edgeR (modo
  `legacy`) usa o maior entre o p do teste F e o de um teste de Poisson. O nosso não faz isso, e
  já dizia que não. São os 2 clusters de 7.600 que ainda diferem do `ql_matched` (p 0,028 contra
  0,110 num, 0,55 contra 0,80 no outro, com a mesma estatística F). O edgeR 4 no modo padrão
  também não aplica esse limite.
- *Clusters sem nenhuma célula.* O nosso não testa; o edgeR devolve p = 1 e conta o cluster nos
  graus de liberdade somados. Na tabela do dia 7 isso deixa 88 graus no nosso e 96 no edgeR, e é
  a diferença de 0,57 dex que resta (p 9e-20 contra 2,5e-20).
- *Um cluster degenerado.* Na tabela dos três dias há um cluster com 1.812 células em uma
  amostra e zero nas outras 29. O log2FC é 5,3 no nosso e 4,1 no `lrt_matched`, e a dispersão
  5,7 contra 7,4: a nossa grade de dispersões termina em 8 e a interpolação é outra. Nenhum dos
  dois números significa muita coisa.
- *Contra os padrões do edgeR 4 e do diffcyt* a diferença é de método, não de conta. O
  `estimateDisp` estima os graus de liberdade da priori das dispersões (mediana 6,3 nas tabelas
  simuladas, de 2,5 a infinito), onde o nosso LRT os fixa em 10 como o
  `estimateGLMTagwiseDisp`. O `glmQLFit` do edgeR 4 usa dispersões com tendência e deviances
  ajustadas, e o nosso QL usa a dispersão comum à maneira anterior à versão 4. Daí 15% a 22% dos
  clusters com p diferente em mais de 0,1 dex, e uma descoberta em dez (LRT contra diffcyt) a
  uma em cinco (QL contra o padrão do edgeR 4) que não coincide. A maior diferença nas simuladas
  (21 dex: p 1,4e-10 no nosso, 5e-32 no diffcyt) é uma tabela em que o `estimateDisp` estimou
  1,3 grau de liberdade para a priori e deu ao cluster uma dispersão de 0,019, onde a priori
  fixa em 10 dá 0,069. Nas tabelas reais as três
  descobertas que só o nosso LRT fez são clusters raros ou ausentes de um grupo (p 0,002 contra
  0,036; 0,006 contra 0,050; 0,001 contra 0,034), justamente onde a dispersão decide.

**Uma observação sobre o edgeR:** `estimateDisp(trend.method = "none")` para com erro
("non-conformable arrays") quando alguma linha tem menos de 5 contagens no total. Nessas tabelas
a linha foi retirada antes de chamar o edgeR e o diffcyt, e não foi comparada. O filtro padrão
do diffcyt (3 células em metade das amostras) já a teria retirado.

---

## 2. PhenoGraph (`milc/graphcluster.py`) e o pacote dos autores

### O que foi comparado

Os mesmos eventos, com k = 30 e resolução 1, três corridas de cada programa. O pacote foi rodado
das duas maneiras que oferece: com o seu Louvain (o programa em C++ de Blondel et al. que vem
dentro dele, que fica com a melhor de muitas partidas aleatórias e não aceita semente) e com
Leiden (leidenalg). O nosso foi chamado só por `graphcluster.phenograph(X, k=, resolution=,
seed=)`. **No commit em que isto foi medido o nosso usa Louvain**; o script pode ser rodado de
novo depois da troca por Leiden sem mudar nada. A rodada com Leiden, só nas misturas, está no
fim desta seção.

Medidas: número de clusters; ARI e NMI entre as partições dos dois programas (média sobre 3 × 3
corridas); ARI entre corridas do mesmo programa, que é a régua; modularidade. Como cada programa
monta o seu grafo, a modularidade das partições de todos foi recalculada no mesmo grafo, o do
pacote, com o mesmo código.

### Dados

- Misturas de 6 gaussianas em 8 dimensões (40%, 25%, 18%, 10%, 6% e 1% dos eventos), as de
  `tests/test_graphcluster.py`, com os centros a 10, 6 e 4,5 desvios-padrão; 6.000 eventos, e
  20.000 no caso de 10 desvios.
- 20.000 eventos CD45+ reais, sorteados das 30 amostras, seis marcadores em escala biexponencial.

### Concordância medida

| Dados | Programa | Clusters (3 corridas) | Modularidade no grafo do pacote | ARI entre as próprias corridas | ARI com a verdade | Segundos por corrida |
|---|---|---|---|---|---|---|
| Mistura, 10 DP, 6 mil | nosso | 6 / 6 / 6 | 0,7616 | 1,00 | 1,00 | 3,1 |
| | pacote, Louvain | 7 / 7 / 7 | 0,7637 | 0,87 | 0,79 | 39 |
| | pacote, Leiden | 6 / 6 / 6 | 0,7616 | 1,00 | 1,00 | 27 |
| Mistura, 6 DP, 6 mil | nosso | 6 / 6 / 6 | 0,7499 | 1,00 | 0,99 | 0,8 |
| | pacote, Louvain | 7 / 7 / 7 | 0,7542 | 0,90 | 0,77 | 46 |
| | pacote, Leiden | 6 / 6 / 6 | 0,7500 | 1,00 | 0,99 | 30 |
| Mistura, 4,5 DP, 6 mil | nosso | 7 / 7 / 6 | 0,7127 | 0,79 | 0,69 | 0,9 |
| | pacote, Louvain | 7 / 7 / 7 | 0,7144 | 0,78 | 0,68 | 33 |
| | pacote, Leiden | 7 / 7 / 7 | 0,7205 | 0,77 | 0,66 | 30 |
| Mistura, 10 DP, 20 mil | nosso | 7 / 7 / 7 | 0,7597 | 0,77 | 0,80 | 6,6 |
| | pacote, Louvain | 7 / 8 / 7 | 0,7586 | 0,76 | 0,76 | 183 |
| | pacote, Leiden | 7 / 8 / 6 | 0,7624 | 0,75 | 0,83 | 105 |
| CD45+ reais, 20 mil | nosso | 19 / 18 / 19 | 0,8535 | 0,71 | — | 2,8 |
| | pacote, Louvain | 20 / 20 / 19 | 0,8551 | 0,73 | — | 133 |
| | pacote, Leiden | 20 / 20 / 18 | 0,8567 | 0,64 | — | 135 |

Um programa contra o outro:

| Dados | Nosso × pacote Louvain (ARI / NMI) | Nosso × pacote Leiden | Pacote Louvain × pacote Leiden |
|---|---|---|---|
| Mistura, 10 DP, 6 mil | 0,79 / 0,92 | 1,00 / 1,00 | 0,79 / 0,92 |
| Mistura, 6 DP, 6 mil | 0,78 / 0,91 | 1,00 / 0,99 | 0,78 / 0,91 |
| Mistura, 4,5 DP, 6 mil | 0,76 / 0,81 | 0,75 / 0,80 | 0,75 / 0,80 |
| Mistura, 10 DP, 20 mil | 0,77 / 0,85 | 0,78 / 0,87 | 0,76 / 0,86 |
| CD45+ reais, 20 mil | 0,71 / 0,81 | 0,64 / 0,78 | 0,66 / 0,78 |

Nenhum programa deixou eventos sem cluster (o pacote descarta clusters com menos de 10 eventos;
não houve nenhum).

### Onde diferem e por quê

- **A concordância entre os programas é a mesma que a de um programa consigo mesmo.** Nos
  eventos reais o nosso concorda com o Louvain do pacote com ARI 0,71; duas corridas do nosso
  concordam com 0,71 e duas do pacote com 0,73. Otimizar modularidade tem muitas soluções quase
  iguais em valor e diferentes no corte das populações grandes; nenhum dos programas é
  reprodutível além disso.
- **O Louvain do pacote divide uma população a mais nas misturas bem separadas** (7 clusters
  onde há 6, cortando a maior em duas). Não é erro: essa partição tem modularidade maior (0,7637
  contra 0,7616), e o pacote a acha porque fica com a melhor de dezenas de partidas, onde o
  nosso faz uma. É o limite de resolução da modularidade que o nosso código já documenta, e
  aparece no nosso também com 20 mil eventos (7 clusters). Com a verdade conhecida, isso dá ao
  nosso ARI maior nesses casos (1,00 contra 0,79), por parar antes, não por ser melhor
  otimizador.
- **A modularidade do nosso fica de 0,002 a 0,004 abaixo da do Louvain do pacote** no mesmo grafo em quatro
  dos cinco casos, pela mesma razão: uma partida contra a melhor de muitas. No caso de 20 mil
  eventos ficou 0,001 acima do Louvain do pacote.
- **Os grafos não são o mesmo grafo.** No nosso, cada evento conta na própria vizinhança e a
  aresta fica com o maior dos dois pesos; no pacote o evento não conta e a aresta fica com a
  média. Por isso a modularidade que o nosso relata (no seu grafo) é de 0,004 a 0,012
  menor do que a da mesma partição no grafo do pacote. As partições, como se vê, passam de um
  grafo para o outro.
- **A modularidade que o pacote relata com Leiden não é comparável** (0,66 a 0,80 onde a mesma
  partição, recalculada, dá 0,72 a 0,86): o leidenalg normaliza de outra forma o grafo que o
  pacote lhe passa. Por isso a coluna recalculada.
- **Tempo.** O nosso levou de 1 a 7 segundos por corrida e o pacote de 27 a 183, nesta máquina
  e com outros processos rodando ao mesmo tempo. Parte é o pacote iniciar processos para o
  Jaccard e, no Louvain, rodar um executável Intel sob Rosetta. Não é uma medida de desempenho
  cuidadosa.

Nada aqui indica erro no nosso código. A troca de Louvain por Leiden, em curso, deve aproximar
o nosso da coluna "pacote, Leiden"; rode o script de novo para ver.

---

### Depois da troca por Leiden

O nosso passou a usar Leiden por padrão no mesmo dia. O script foi rodado de novo, sem mudança,
**só nas misturas simuladas**; os 20.000 eventos reais não foram repetidos.

| Dados | Programa | Clusters (3 corridas) | Modularidade no grafo do pacote | ARI entre as próprias corridas | ARI com a verdade | Segundos por corrida |
|---|---|---|---|---|---|---|
| Mistura, 10 DP, 6 mil | nosso (Leiden) | 6 / 6 / 6 | 0,7616 | 1,00 | 1,00 | 0,8 |
| | pacote, Louvain | 7 / 7 / 7 | 0,7636 | 0,88 | 0,79 | 30 |
| | pacote, Leiden | 6 / 6 / 6 | 0,7616 | 1,00 | 1,00 | 25 |
| Mistura, 6 DP, 6 mil | nosso (Leiden) | 6 / 6 / 6 | 0,7500 | 1,00 | 0,99 | 0,6 |
| | pacote, Louvain | 7 / 7 / 7 | 0,7534 | 0,83 | 0,78 | 35 |
| | pacote, Leiden | 6 / 6 / 6 | 0,7500 | 1,00 | 0,99 | 27 |
| Mistura, 4,5 DP, 6 mil | nosso (Leiden) | 7 / 7 / 6 | 0,7171 | 0,76 | 0,76 | 0,7 |
| | pacote, Louvain | 7 / 7 / 7 | 0,7146 | 0,74 | 0,68 | 33 |
| | pacote, Leiden | 7 / 7 / 7 | 0,7205 | 0,77 | 0,66 | 27 |
| Mistura, 10 DP, 20 mil | nosso (Leiden) | 7 / 7 / 7 | 0,7664 | 0,79 | 0,79 | 4,0 |
| | pacote, Louvain | 7 / 7 / 7 | 0,7590 | 0,76 | 0,79 | 167 |
| | pacote, Leiden | 7 / 8 / 6 | 0,7624 | 0,75 | 0,83 | 80 |

O nosso Leiden contra o Leiden do pacote: ARI 1,00 a 10 DP e 0,997 a 6 DP com 6 mil eventos;
0,75 a 4,5 DP e 0,78 com 20 mil eventos, onde o próprio pacote concorda consigo mesmo em 0,75 a
0,77. Com 20 mil eventos o nosso passou a dividir a maior população (7 clusters, como o pacote),
o que o Louvain anterior não fazia: um otimizador melhor acha a partição de maior modularidade,
que é a dividida. Nos 3.488 eventos reais de teste, medidos à parte, o nosso e o pacote dão 15
clusters cada, com ARI de 0,82 a 0,89 entre eles.

## 3. Controle de qualidade da aquisição (`milc/qc.py`), o PeacoQC e o flowAI

Esta seção foi medida duas vezes no mesmo dia. A primeira medida achou que o nosso deixava passar
desvios moderados e sustentados; `qc.check` ganhou então uma segunda passada, lenta, para eles
(parâmetro `drift`), e tudo foi medido de novo. As tabelas trazem as duas medidas, "antes" e
"depois". Os números do PeacoQC e do flowAI são os mesmos nas duas.

### O que foi comparado

Quais eventos cada programa remove de uma mesma aquisição.

- **Nosso**: `qc.check` chamado como o programa chama (`Session.quality`): todos os parâmetros
  menos Time, fluorescência compensada, cada parâmetro na escala em que é exibido, padrões
  (`mad = 6`, `consecutive = 5` e, depois, `drift = 0,25`).
- **PeacoQC**: recebe exatamente os mesmos valores que o nosso julgou, na ordem do arquivo, nos
  seus padrões (`MAD = 6`, `IT_limit = 0.6`, `consecutive_bins = 5`). Rodado inteiro (árvore de
  isolamento e depois MAD) e só com o passo MAD, que é o que o nosso imita.
- **flowAI**: recebe os valores como adquiridos, que é o que ele espera, com as checagens de
  vazão e de sinal (`remove_from = "FR_FS"`). Conferido em um arquivo (3Case2) que ele dá o mesmo
  resultado lendo o FCS diretamente (5,994% nos dois casos).

O PeacoQC não tem checagem de vazão: ele só olha o sinal. O flowAI tem as duas.

### Dados

- Cinco aquisições sintéticas como as de `tests/test_qc.py` (60 s a 1.000 eventos/s, dois
  parâmetros de dispersão e dois de fluorescência com duas populações), uma limpa e quatro com
  um defeito cada.
- Os 32 arquivos reais, inteiros (de 13,8 mil a 247 mil eventos, 12 parâmetros).
- Um arquivo real limpo (14Case5: o nosso remove 0,4%, antes e depois, e o PeacoQC 0,8%) com um
  defeito escrito nos valores brutos, três vezes: um degrau (fluorescência × 0,5 em 40–50% dos
  eventos), uma deriva (fluorescência caindo até × 0,5 ao longo dos últimos 30%) e um
  entupimento (os eventos de 60–63% chegando dez vezes mais devagar, sem mudar o sinal).

### Onde o defeito é conhecido

Percentual dos eventos do defeito que cada programa remove, e percentual dos outros eventos que
ele remove junto:

| Aquisição | Defeito (% dos eventos) | Nosso antes: do defeito | Antes: do resto | Nosso depois: do defeito | Depois: do resto | PeacoQC: do defeito | PeacoQC: do resto | PeacoQC só MAD: do defeito | Só MAD: do resto | flowAI: do defeito | flowAI: do resto |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Sintética, entupimento (vazão a 3% por 10 s) | 0,56 | 100 | 0,34 | 100 | 0,34 | 0 | 0 | 0 | 0 | 4,6 | 31,5 |
| Sintética, degrau (PE × 1,8 nos últimos 15 s) | 24,9 | 100 | 0,16 | 100 | 0,16 | 100 | 1,16 | 100 | 2,28 | 28,9 | 11,1 |
| Sintética, deriva (FITC até × 3 nos últimos 20 s) | 33,5 | 81,1 | 0 | 97,3 | 0 | 97,6 | 0 | 95,1 | 0 | 24,6 | 11,5 |
| Sintética, queda (tudo × 0,5 por 2 s) | 3,5 | 98,9 | 0,05 | 98,9 | 0,05 | 100 | 1,99 | 100 | 5,86 | 4,1 | 14,3 |
| 14Case5 + degrau | 10,0 | 99,2 | 0,55 | 100 | 0,69 | 99,2 | 1,36 | 100 | 3,52 | 3,7 | 7,8 |
| 14Case5 + deriva | 30,0 | 53,5 | 0,57 | 76,2 | 0,57 | 74,2 | 0,58 | 74,2 | 0,58 | 82,3 | 7,2 |
| 14Case5 + entupimento | 3,0 | 100 | 0,63 | 100 | 0,63 | 0 | 0,84 | 0 | 0,84 | 8,6 | 12,6 |

Na aquisição sintética limpa o nosso não remove nada, antes e depois; o PeacoQC remove 0,6% e o
flowAI 13,2%.

Sobreposição dos conjuntos removidos (Jaccard) entre o nosso e o PeacoQC, antes → depois: 0,97 →
0,97 no degrau sintético, 0,92 → 0,90 no degrau injetado, 0,83 → 1,00 na deriva sintética, 0,72
→ 0,98 na injetada, 0,65 na queda, 0 e 0,10 nos entupimentos.

### Arquivos reais

Os 32 juntos:

| Programa | Medida | Nosso remove (% de todos os eventos) | Ele remove | Mediana por arquivo: nosso | Mediana: ele | Jaccard sobre todos os eventos | Do que o nosso remove, ele também remove (%) | Do que ele remove, o nosso também (%) | Arquivos: os dois abaixo de 1% | Os dois com 1% ou mais | Só o nosso com 1% ou mais | Só ele |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PeacoQC | antes | 1,37 | 5,13 | 0,40 | 2,48 | 0,19 | 75,3 | 20,1 | 8 | 10 | 3 | 11 |
| PeacoQC | depois | 5,25 | 5,13 | 0,52 | 2,48 | 0,38 | 54,9 | 56,1 | 8 | 11 | 3 | 10 |
| PeacoQC só MAD | antes | 1,37 | 4,60 | 0,40 | 2,09 | 0,21 | 75,3 | 22,5 | 8 | 10 | 3 | 11 |
| PeacoQC só MAD | depois | 5,25 | 4,60 | 0,52 | 2,09 | 0,41 | 54,9 | 62,6 | 8 | 11 | 3 | 10 |
| flowAI | antes | 1,37 | 11,97 | 0,40 | 7,02 | 0,08 | 73,6 | 8,4 | 0 | 13 | 0 | 19 |
| flowAI | depois | 5,25 | 11,97 | 0,52 | 7,02 | 0,17 | 46,8 | 20,5 | 0 | 14 | 0 | 18 |

Arquivo por arquivo (percentual dos eventos removido). A passada lenta mudou seis arquivos; os
outros 26 têm exatamente os mesmos eventos removidos antes e depois.

| Arquivo | Eventos | Nosso antes | Nosso depois | Depois: vazão | Depois: sinal | Depois: deriva | PeacoQC | PeacoQC só MAD | flowAI | Jaccard nosso/PeacoQC antes | Depois | Jaccard nosso/flowAI depois |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 14Case1 | 166.402 | 6,9 | 22,8 | 0,0 | 6,1 | 22,0 | 14,1 | 14,1 | 69,0 | 0,44 | 0,62 | 0,33 |
| 14Case2 | 119.308 | 3,0 | 4,2 | 0,2 | 3,0 | 4,2 | 3,1 | 3,1 | 7,8 | 0,96 | 0,75 | 0,20 |
| 14Case3 | 122.291 | 0,2 | 16,5 | 0,2 | 0,0 | 16,2 | 3,5 | 3,5 | 13,9 | 0,07 | 0,20 | 0,03 |
| 14Case4 | 140.573 | 1,0 | 1,0 | 0,8 | 1,0 | 0,0 | 1,5 | 1,5 | 11,0 | 0,68 | 0,68 | 0,08 |
| 14Case4` | 109.370 | 2,1 | 11,3 | 1,5 | 0,0 | 9,4 | 0,5 | 0,5 | 14,5 | 0,09 | 0,04 | 0,18 |
| 14Case5 | 247.019 | 0,4 | 0,4 | 0,4 | 0,0 | 0,0 | 0,8 | 0,8 | 7,3 | 0,49 | 0,49 | 0,05 |
| 14Ctrl1 | 240.680 | 1,1 | 7,1 | 0,5 | 0,8 | 6,6 | 13,8 | 13,8 | 6,6 | 0,06 | 0,49 | 0,21 |
| 14Ctrl2 | 144.357 | 0,4 | 0,4 | 0,2 | 0,2 | 0,0 | 1,7 | 1,7 | 6,0 | 0,10 | 0,10 | 0,03 |
| 14Ctrl3 | 217.920 | 0,4 | 0,4 | 0,4 | 0,0 | 0,0 | 0,0 | 0,0 | 7,3 | 0,00 | 0,00 | 0,05 |
| 14Ctrl4 | 194.677 | 0,2 | 0,2 | 0,2 | 0,0 | 0,0 | 1,6 | 1,6 | 6,6 | 0,12 | 0,12 | 0,01 |
| 14Ctrl5 | 13.816 | 7,7 | 7,7 | 2,3 | 4,3 | 0,0 | 2,3 | 2,3 | 4,9 | 0,30 | 0,30 | 0,03 |
| 3Case1 | 26.242 | 4,6 | 4,6 | 0,0 | 4,6 | 0,0 | 3,8 | 3,8 | 14,3 | 0,83 | 0,83 | 0,24 |
| 3Case2 | 44.511 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 6,0 | 1,00 | 1,00 | 0,00 |
| 3Case3 | 42.226 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 6,1 | 1,00 | 1,00 | 0,00 |
| 3Case4 | 39.048 | 0,8 | 0,8 | 0,4 | 0,8 | 0,0 | 0,1 | 0,1 | 7,5 | 0,00 | 0,00 | 0,10 |
| 3Case5 | 48.068 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 4,7 | 1,00 | 1,00 | 0,00 |
| 3Ctrl1 | 16.337 | 2,7 | 2,7 | 2,7 | 2,7 | 0,0 | 5,1 | 5,1 | 9,0 | 0,52 | 0,52 | 0,30 |
| 3Ctrl2 | 24.263 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 16,5 | 16,5 | 5,3 | 0,00 | 0,00 | 0,00 |
| 3Ctrl3 | 30.623 | 1,5 | 1,5 | 0,6 | 1,0 | 0,0 | 2,9 | 2,9 | 5,2 | 0,54 | 0,54 | 0,18 |
| 3Ctrl4 | 20.043 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 4,0 | 4,0 | 5,1 | 0,00 | 0,00 | 0,00 |
| 3Ctrl5 | 16.765 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 1,6 | 1,6 | 6,8 | 0,00 | 0,00 | 0,00 |
| 7Case1 | 47.784 | 0,6 | 0,6 | 0,6 | 0,0 | 0,0 | 12,1 | 10,0 | 5,2 | 0,03 | 0,03 | 0,05 |
| 7Case2 | 53.920 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 3,7 | 6,5 | 8,2 | 0,00 | 0,00 | 0,00 |
| 7Case3 | 16.139 | 2,4 | 2,4 | 2,4 | 0,0 | 0,0 | 0,0 | 0,0 | 20,7 | 0,00 | 0,00 | 0,01 |
| 7Case3` | 65.720 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 2,7 | 1,9 | 6,9 | 0,00 | 0,00 | 0,00 |
| 7Case4 | 95.972 | 0,4 | 0,4 | 0,2 | 0,4 | 0,0 | 0,5 | 0,5 | 8,7 | 0,00 | 0,00 | 0,00 |
| 7Case5 | 130.097 | 4,5 | 4,5 | 2,6 | 3,5 | 2,9 | 3,5 | 3,5 | 8,9 | 0,78 | 0,78 | 0,33 |
| 7Ctrl1 | 27.207 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 8,1 | 8,1 | 29,5 | 0,00 | 0,00 | 0,00 |
| 7Ctrl2 | 44.866 | 0,3 | 0,3 | 0,3 | 0,0 | 0,0 | 0,8 | 0,8 | 6,2 | 0,00 | 0,00 | 0,00 |
| 7Ctrl3 | 40.307 | 2,6 | 2,6 | 0,4 | 2,2 | 0,0 | 0,0 | 0,0 | 5,1 | 0,00 | 0,00 | 0,08 |
| 7Ctrl4 | 60.261 | 0,0 | 0,0 | 0,0 | 0,0 | 0,0 | 31,1 | 7,1 | 7,1 | 0,00 | 0,00 | 0,00 |
| 7Ctrl5 | 103.643 | 1,0 | 32,7 | 0,8 | 0,2 | 32,7 | 19,4 | 19,4 | 6,8 | 0,05 | 0,58 | 0,13 |

(Jaccard 1,00 onde nenhum dos dois remove nada. Um bin achado por mais de uma checagem do nosso
conta em cada coluna.)

### O que mudou em `qc.py`, e por que assim

A primeira medida mostrou a falha: o nosso julgava a mediana corrida de 5 bins contra a mediana
do arquivo, com um limite em MADs dessa mesma curva. Um deslocamento de 0,3 do desvio dos eventos
fica dentro do ruído de um bin de 150 a 500 eventos, e a mediana corrida de 5 bins não o tira de
lá. Em 7Ctrl5 o primeiro terço do arquivo tem quatro canais deslocados de 0,25 a 0,35 do desvio
e o nosso removia 1,0%.

O que foi posto (`drift`, padrão 0,25; 0 desliga):

1. As mesmas medianas de bin passam por uma mediana corrida longa: a janela tem 10.000 eventos,
   e nunca mais que um quarto dos bins.
2. O nível de referência é o do trecho estável mais longo: o meio da faixa de valores mais
   estreita que contém metade dos bins. Não é a mediana.
3. Um trecho é deriva quando a curva lenta se afasta desse nível mais que `mad` MADs da própria
   curva lenta *e* mais que `drift` vezes o desvio dos eventos do parâmetro.
4. Do trecho achado remove-se o afastamento inteiro, desde onde a curva lenta está a 0,7 do
   limite, e não só a parte além dele.

Cada uma dessas escolhas foi medida contra a alternativa, nos mesmos 40 casos (percentual dos
eventos removido; "reais" é o total dos 32 arquivos):

| Variante | Reais | 7Ctrl5 | 14Case1 | 3Ctrl2 | 7Case1 | 7Ctrl4 | Deriva sintética (do defeito) | Deriva injetada (do defeito) |
|---|---|---|---|---|---|---|---|---|
| Antes (sem passada lenta) | 1,37 | 1,0 | 6,9 | 0 | 0,6 | 0 | 81,1 | 53,5 |
| **Como ficou** | 5,25 | 32,7 | 22,8 | 0 | 0,6 | 0 | 97,3 | 76,2 |
| Nível pela mediana do arquivo | 3,85 | 18,5 | 22,2 | 0 | 0,6 | 0 | 95,1 | 75,5 |
| Sem o piso (`drift = 0,001`) | 10,44 | 39,6 | 22,8 | 25,8 | 24,7 | 17,2 | 100 | 86,9 |
| Sem o MAD: limite fixo de 0,25 do desvio | 24,17 | 32,7 | 100 | 0 | 0,6 | 0 | 97,3 | 100 |
| Sem o MAD: limite fixo de 0,4 do desvio | 12,30 | 19,5 | 80,6 | 0 | 0,6 | 0 | 94,4 | 66,1 |
| Só a parte além do limite | 4,04 | 20,5 | 20,8 | 0 | 0,6 | 0 | 94,4 | 65,5 |
| Desde 0,5 do limite | 5,95 | 33,3 | 26,2 | 0 | 0,6 | 0 | 100 | 87,5 |
| Janela de 2.500 eventos | 4,11 | 17,5 | 22,8 | 0 | 0,6 | 0 | 94,4 | 61,5 |
| Janela de 5.000 eventos | 4,74 | 24,3 | 24,2 | 0 | 0,6 | 0 | 95,1 | 73,5 |
| Janela de 20.000 eventos | 5,94 | 26,9 | 25,0 | 0 | 0,6 | 0 | 96,6 | 76,2 |
| `drift = 0,15` | 6,31 | 32,7 | 22,8 | 0 | 24,4 | 0 | 100 | 78,2 |
| `drift = 0,4` | 3,40 | 19,5 | 22,6 | 0 | 0,6 | 0 | 94,4 | 62,8 |
| `mad = 5` (as duas passadas) | 6,14 | 33,0 | 27,0 | 1,2 | 0,9 | 0,2 | 97,3 | 76,2 |
| `mad = 8` (as duas passadas) | 2,99 | 21,7 | 20,2 | 0 | 0 | 0 | 96,6 | 73,5 |
| PeacoQC | 5,13 | 19,4 | 14,1 | 16,5 | 12,1 | 31,1 | 97,6 | 74,2 |

O que cada linha diz:

- *Nível pela mediana*: a mediana é puxada para o trecho deslocado. De 7Ctrl5 saem 18,5% em vez
  de 32,7%, e o ressalto de 14Case3 (16,5%) não é visto.
- *Sem o piso* (`drift` quase zero): é o critério do PeacoQC, MADs de uma curva alisada e mais
  nada. Um arquivo cuja curva quase não se mexe tem um MAD minúsculo e perde os bins onde ela
  mais se mexeu: 3Ctrl2, 7Case1 e 7Ctrl4 perdem 25,8%, 24,7% e 17,2%, como no PeacoQC, sem que nenhuma
  mediana de canal se desloque mais de 0,05 do desvio. O piso é o que impede isso.
- *Sem o MAD* (só o piso, como limite fixo): as aquisições reais oscilam mais que 0,25 do desvio
  ao longo do arquivo todo, e um limite fixo corta as que oscilam: 24% de todos os eventos
  reais, 14Case1 e 14Case3 inteiros, e um terço do resto do arquivo com o degrau injetado. Com
  0,4, ainda 12%. O MAD da curva lenta é o limite do próprio arquivo.
- *Só a parte além do limite* (sem o item 4): de uma rampa só sai o trecho final, e de 7Ctrl5
  saem 20,5% quando o deslocamento dura um terço do arquivo.
- *Desde 0,5 do limite*: pega mais da deriva injetada, e começa a remover antes do defeito (na
  deriva sintética, 1,3% do resto; com 0,7, nada).
- *Janela*: menor, vê menos das derivas; maior, tira mais dos arquivos reais sem achar mais do
  defeito conhecido (com 20.000, 34% de 14Case3).
- *`drift` e `mad`*: são os dois limites que estão no diálogo. Com `drift = 0,15` 7Case1 perde
  24,4%: o piso de 0,25 tem folga para os outros 31 arquivos, não para esse.

Uma variante foi medida só no protótipo, antes do item 4: o desvio dos eventos *da população*,
e não do parâmetro, como unidade do piso. As populações empilhadas no topo da escala têm desvio
quase nulo, qualquer movimento delas passa do piso, e 3Ctrl2 perdia 8,7%.

Em simulação, com as aquisições de `tests/test_qc.py` (percentual do defeito removido / do
resto, média de 10 sementes; 60 s a 1.000 eventos/s):

| Defeito | Antes | Depois |
|---|---|---|
| FSC mais alto nos primeiros 20%, em 0,20 do desvio | 0 / 0 | 0 / 0 |
| idem, 0,25 | 0 / 0 | 56,0 / 0 |
| idem, 0,30 | 0 / 0 | 97,9 / 0,08 |
| idem, 0,40 | 8,8 / 0 | 99,8 / 0,40 |
| idem, 0,50 | 45,8 / 0 | 100 / 0,45 |
| idem, 0,75 | 99,6 / 0,02 | 100 / 0,42 |
| Só a população mais brilhante de FITC (40% das células) × 1,10 nos primeiros 20% | 0 / 0 | 95,5 / 0,39 |
| idem, × 1,15 | 0,1 / 0 | 99,9 / 1,58 |
| idem, × 1,30 | 92,5 / 0,01 | 100 / 1,92 |
| Rampa em FSC nos últimos 30%, até 0,5 do desvio | 3,0 / 0 | 64,6 / 0 |
| idem, até 1 desvio | 42,6 / 0 | 82,2 / 0 |
| idem, até 2 desvios | 68,7 / 0 | 90,3 / 0 |
| FSC mais alto em 0,3 do desvio durante 5% do arquivo | 13,5 / 0 | 13,5 / 0 |
| idem, 0,5 | 92,5 / 0,02 | 92,5 / 0,02 |

Em 520 aquisições limpas simuladas (de 15 a 250 s, de 300 a 1.000 eventos/s, de 10% a 48% das
células na população brilhante) a passada lenta não removeu nenhum evento: as 6 que perdem algo
(0,4 a 0,75 por mil dos eventos do seu grupo) perdem o mesmo antes e depois. Uma senoide em FSC
ao longo do arquivo todo (1,5 período, até 0,8 do desvio) não tira nada.

### Onde diferem e por quê

1. **Desvios lentos e moderados: era a discrepância principal, e fechou.** Nos arquivos reais o
   nosso remove agora 5,25% dos eventos (antes 1,37%) e o PeacoQC 5,13%; mas não são os mesmos
   eventos: do que cada um remove, o outro remove pouco mais da metade (55% e 56%; antes, 75% e
   20%), pelos itens 2 e 4. Os dois casos abertos na primeira medida:
   - *7Ctrl5*: a mediana de PerCP-Cy5.5 está deslocada em 0,35 do desvio dos eventos, e as de
     PE, APC e BV421 em 0,25 a 0,30, no começo do arquivo (a primeira medida olhou os 17% que o
     PeacoQC remove; a curva lenta mostra que vai até um terço). O PeacoQC remove os primeiros
     17% (19,4% ao todo); o nosso removia 1,0% e remove agora os primeiros 31% do tempo (32,7%
     dos eventos), que é até onde o deslocamento vai (o degrau está em 33% dos bins).
   - *14Case1*: o PeacoQC remove 14,1%; o nosso removia 6,9% e remove agora 22,8%, de 57% a 95%
     do tempo: o trecho em que APC, FSC e SSC afundam.

   Nos eventos que só a passada lenta remove, as medianas estão deslocadas, em relação aos
   eventos mantidos, nesta fração do desvio dos eventos do canal:

   | Arquivo | Nosso antes → depois (%) | Canais mais deslocados |
   |---|---|---|
   | 14Case1 | 6,9 → 22,8 | APC-Cy7 1,23; PE-Cy7 0,69; APC 0,59; SSC-H 0,43 |
   | 14Case2 | 3,0 → 4,2 | PE 0,36; FITC 0,26 |
   | 14Case3 | 0,2 → 16,5 | FITC 0,29; SSC-H 0,29; SSC-A 0,27; BV421 0,22 |
   | 14Case4` | 2,1 → 11,3 | SSC-A 0,44; PE 0,33; PerCP-Cy5.5 0,31; SSC-H 0,31 |
   | 14Ctrl1 | 1,1 → 7,1 | FITC 0,27; PE 0,17; SSC-H 0,16 |
   | 7Ctrl5 | 1,0 → 32,7 | PerCP-Cy5.5 0,30; BV421 0,30; PE 0,29; BV510 0,26 |

   A causa provável apontada na primeira medida era o alisamento (o PeacoQC calcula a mediana e o
   MAD de uma curva alisada por `smooth.spline(spar = 0.5)`) e o tamanho dos bins (500 ou 1.000
   eventos, com sobreposição de metade, contra 150 a 495 nos nossos). O alisamento foi testado e
   era isso: uma janela longa com o MAD da curva alisada fecha a diferença. Os bins não foram
   aumentados. O comentário em `qc.py` que dizia que 150 e 500 eram "os limites do PeacoQC" foi
   corrigido: 150 é só o piso dele, nunca atingido, e o bin dele tem o dobro de 1/500 do arquivo
   arredondado para cima em múltiplos de 500.
2. **Remoções grandes do PeacoQC sem deslocamento visível das medianas: o nosso continua sem
   fazê-las.** Em 3Ctrl2 (16,5%), 7Case1 (12,1%) e 7Ctrl4 (31,1%) o nosso remove 0, 0,6% e 0,
   antes e depois. Nos eventos que o PeacoQC remove, nenhuma mediana de canal difere do resto em
   mais de 0,05 do desvio dos eventos (menos de 0,6 da dispersão das medianas dos bins). O
   próprio PeacoQC atribui essas remoções à posição de picos de densidade (o segundo pico de PE
   em 7Ctrl4; o pico de PE-Cy7 e o segundo de BV421 em 3Ctrl2) e, em 7Ctrl4 e 7Case1, à árvore de
   isolamento. O PeacoQC segue a moda de cada pico e o nosso a mediana de cada população. O
   nosso não segue populações com menos de 20 eventos por bin e não tem árvore de isolamento.
   Sem a verdade não dá para dizer quem tem razão; a medida diz só que o canal inteiro não se
   moveu. A passada lenta sem o piso faz o mesmo que o PeacoQC nesses três arquivos (tabela
   acima); com o piso, não.
3. **Vazão.** O PeacoQC não a examina: não removeu nada dos dois entupimentos, que o nosso
   removeu inteiros. O flowAI examina, de modo muito mais estrito que o nosso: marcou de 4,7% a
   9% dos eventos em 25 dos 32 arquivos reais e mais nos outros (até 69%). Em 3Case2, o único
   aberto, os 5,994% vêm todos da checagem de vazão, em dezenas de trechos curtos. E não acertou os entupimentos (4,6% e 8,6% do defeito, removendo
   31% e 13% do resto). O nosso só reage a uma vazão fora de um fator 1,5 da habitual ou a 6
   MADs.
4. **O nosso remove e o PeacoQC não, ou remove menos.** Já era assim em 14Ctrl5 (7,7% contra
   2,3%), 7Ctrl3 (2,6% contra 0) e, pela vazão, 7Case3 (2,4%) e 14Case4` (2,1% contra 0,5%). Com a
   passada lenta: 14Case3 (16,5% contra 3,5%: um ressalto de SSC de 41% a 58% do tempo, do qual o
   PeacoQC remove o meio), 14Case4` (11,3% contra 0,5%: PE e PerCP-Cy5.5 mais altos nos
   primeiros 8%), e 7Ctrl5 e 14Case1, onde o nosso remove o afastamento inteiro e o PeacoQC a
   parte mais afastada. Em 14Ctrl1 é o contrário: 7,1% contra 13,8% (o PeacoQC remove os
   primeiros 11%, o nosso os primeiros 7%).
5. **Falsos alarmes.** Onde o defeito é conhecido, o nosso remove menos do que não devia: de 0 a
   0,7% do resto (antes, de 0 a 0,6%), contra 0 a 2,0% do PeacoQC inteiro e até 5,9% do passo
   MAD sozinho.
6. **O flowAI quase não viu os degraus e a queda** (3,7% a 29% do defeito) e removeu de 8% a 14%
   do resto. Ele trabalha nos valores brutos, sem transformação; é provável que isso explique o
   degrau na fluorescência, mas não foi investigado. Para o flowAI esta comparação é só um
   retrato do uso padrão.

### O que piorou, e o que continua sem ser visto

- No degrau injetado o nosso remove 0,69% do resto, onde removia 0,55%.
- Dois arquivos reais perdem mais de um quinto (7Ctrl5 32,7%, 14Case1 22,8%) e dois mais de um
  décimo (14Case3, 14Case4`). Em todos há deslocamento medido (tabela do item 1), mas não há a
  verdade: ninguém sabe se esses eventos estavam errados. O diálogo avisa quando passa de um
  quinto, e os limites estão nele para mudar.
- Um defeito moderado e curto continua invisível: a passada rápida não tem sensibilidade e a
  lenta não tem resolução. Em simulação, 0,3 do desvio durante 5% do arquivo: 13,5% do
  defeito removido, antes e depois.
- Uma deriva que atravessa o arquivo inteiro sem trecho estável não é removida (em simulação,
  uma senoide de até 0,8 do desvio: nada). É de propósito: não há nível de que se afastar.
- Quando mais da metade do arquivo está deslocada, o que sai é a outra parte (em simulação, 60%
  deslocado: saem os outros 40%).
- O resultado depende da escala em que o parâmetro é julgado. 14Case1 com as escalas do workspace
  do FlowJo dá os 22,8% da tabela; aberto sem o workspace, nas escalas padrão do programa, dá
  13,0%.

Resumo: nos defeitos nítidos o nosso continua igual ou melhor que o PeacoQC e tem a checagem de
vazão que ele não tem; nos desvios sustentados de 0,3 a 0,7 desvio, que eram a falha, está agora
no nível do PeacoQC (97,3% contra 97,6% e 76,2% contra 74,2% do defeito conhecido), sem as
remoções dele em arquivos cujas medianas não se movem.

---

## 4. Normalização de lotes (`milc/normalize.py`) e o CytoNorm

### O que já tinha sido comparado

`tests/test_normalize_r.py` compara as nossas curvas com três linhas de R base (quantis de 0,01
a 0,99, média dos lotes como alvo, `splinefun(method = "monoH.FC")`), que é o miolo do
`QuantileNorm` para um cluster. O pacote CytoNorm em si, com o FlowSOM, não tinha sido usado.

### O que foi comparado agora

O pacote inteiro (CytoNorm 2.0.13, do repositório dos autores), nos mesmos eventos. Os eventos
vão para o R já na escala em que se normaliza, em arquivos FCS de 32 bits que o próprio R
escreve, de modo que os dois programas recebem os mesmos números e nenhuma transformação entra
na comparação.

- **Sem clusters**: o nosso `fit(n_clusters=1)` contra `QuantileNorm.train` /
  `QuantileNorm.normalize`.
- **Com 10 metaclusters**: o nosso `fit(n_clusters=10)` (mapa 10 × 10, o nosso padrão, e também
  15 × 15) contra `CytoNorm.train` / `CytoNorm.normalize` (mapa 15 × 15, o padrão dele, e também
  10 × 10), `nQ = 99`, `scale = FALSE`. Cada programa com duas sementes: os mapas dos dois nunca
  são o mesmo mapa, então a régua é o quanto duas sementes do mesmo programa diferem.

### Dados

- **Simulados, com a verdade conhecida**: três lotes, três populações em quatro marcadores, um
  efeito de lote diferente para cada população (deslocamentos de até 0,6 do desvio-padrão da
  população e escala de 0,9 a 1,15). O lote "a" não tem efeito e é o alvo, de modo que cada
  evento normalizado deveria voltar ao valor que tinha antes do efeito. 20.000 eventos de
  referência por lote; 10.000 a normalizar por lote, com outra composição (20/30/50% em vez de
  50/30/20%).
- **Reais**: o dia como lote (3, 7, 14), os animais Ctrl como referência (12.738, 30.867 e
  68.057 células vivas), os animais Case a normalizar (40.488, 59.910 e 77.923), oito marcadores
  em escala biexponencial (eixo de 0 a 1), alvo a média dos dias.

### Concordância medida

Sem clusters, nos eventos a normalizar, dentro da faixa de 1% a 99% das referências:

| Dados, lote | Correção média feita | Diferença média nosso − CytoNorm | p99 | Máx. | Valores fora da faixa (%) | Diferença média fora | Máx. fora |
|---|---|---|---|---|---|---|---|
| Real, dia 3 | 0,0199 | 0,00001 | 0,0002 | 0,0022 | 2,9 | 0,0067 | 0,092 |
| Real, dia 7 | 0,0166 | 0,00002 | 0,0004 | 0,0021 | 2,8 | 0,0069 | 0,085 |
| Real, dia 14 | 0,0385 | 0,00002 | 0,0005 | 0,0014 | 8,2 | 0,0070 | 0,17 |
| Simulado, b | 0,084 | 0,0003 | 0,0084 | 0,029 | 2,2 | 0,015 | 0,13 |
| Simulado, c | 0,085 | 0,0007 | 0,022 | 0,049 | 2,1 | 0,009 | 0,12 |

(Unidades da escala. No simulado o lote "a" é o alvo e não muda em nenhum dos dois: diferença
zero.)

Com 10 metaclusters, a que distância ficam duas normalizações dos mesmos eventos:

| Par | Simulado: diferença média | Em % da correção média | Real: diferença média | Em % da correção média |
|---|---|---|---|---|
| Nosso, semente 0 × semente 1 | 0,0098 | 18 | 0,0038 | 24 |
| CytoNorm, semente 1 × semente 2 | 0,0093 | 17 | 0,0049 | 31 |
| Nosso (10 × 10) × CytoNorm (15 × 15) | 0,0075 | 14 | 0,0042 | 27 |
| Nosso (10 × 10) × CytoNorm (10 × 10) | 0,0107 | 20 | 0,0061 | 39 |
| Nosso (15 × 15) × CytoNorm (15 × 15) | 0,0082 | 15 | 0,0035 | 23 |
| Nosso com clusters × nosso sem | 0,035 | 65 | 0,019 | 117 |
| CytoNorm com clusters × CytoNorm sem | 0,034 | 61 | 0,019 | 121 |

O que sobra do efeito de lote (`normalize.distance`: distância média entre os quantis dos lotes
e o alvo, em unidades da escala) e, no simulado, o erro contra a verdade:

| Normalização | Simulado: referências | Simulado: erro médio contra o valor verdadeiro | p99 do erro | Real: referências | Real: amostras Case |
|---|---|---|---|---|---|
| Antes | 0,053 | 0,085 | 0,20 | 0,0224 | 0,0300 |
| Nosso, sem clusters | 0,0000 | 0,051 | 0,32 | 0,0000 | 0,0280 |
| CytoNorm, sem clusters | 0,0000 | 0,051 | 0,33 | 0,0000 | 0,0280 |
| Nosso, 10 clusters (10 × 10) | 0,0166 | 0,012 | 0,092 | 0,0135 | 0,0288 |
| Nosso, 10 clusters (15 × 15) | 0,0165 | 0,012 | 0,085 | 0,0134 | 0,0277 |
| CytoNorm, 10 clusters (15 × 15) | 0,0110 | 0,009 | 0,12 | 0,0138 | 0,0273 |
| CytoNorm, 10 clusters (10 × 10) | 0,0233 | 0,022 | 0,096 | 0,0131 | 0,0295 |

(A distância entre as amostras Case não deve ir a zero: elas diferem de um dia para outro por
biologia. Com clusters a distância entre as referências também não vai a zero, porque a medida
olha cada marcador inteiro e a correção é feita dentro dos clusters.)

Até onde um valor foi levado (dados reais, eixo de 0 a 1, 1,43 milhão de valores):

| Normalização | Maior correção de um valor | Menor valor depois | Maior valor depois | Valores levados mais de 10% do eixo para fora da faixa da entrada |
|---|---|---|---|---|
| Nosso, sem clusters | 0,25 | −0,011 | 1,15 | 14 |
| CytoNorm, sem clusters | 0,25 | −0,18 | 1,14 | 231 |
| Nosso, 10 clusters (10 × 10) | 0,17 | −0,040 | 1,04 | 0 |
| Nosso, 10 clusters (15 × 15) | 0,16 | −0,025 | 1,11 | 1 |
| CytoNorm, 10 clusters (15 × 15), semente 1 | 0,99 | −0,15 | 1,51 | 19 |
| CytoNorm, 10 clusters (10 × 10) | 1,12 | −0,89 | 1,55 | 124 |

Na segunda semente do CytoNorm 15 × 15 um valor do dia 7 foi levado a 4,96 (correção de 4,2), que
é a maior diferença entre duas corridas do próprio CytoNorm nesses dados.

### Onde diferem e por quê

- **Dentro da faixa das referências, sem clusters, é a mesma curva.** A diferença média é de um
  milésimo da correção nos dados reais e de 0,3% a 0,8% no simulado. O que resta vem de
  `PchipInterpolator` e do `monoH.FC` do R escolherem as inclinações nos nós de modo um pouco
  diferente, o que aparece onde os nós estão longe um do outro (o vão entre duas populações),
  como o teste antigo já dizia.
- **Fora da faixa, a diferença é de propósito.** Além do primeiro e do último quantil o CytoNorm
  continua a inclinação da ponta da curva, e o nosso continua como deslocamento. Com clusters
  pequenos a ponta é estimada com poucos eventos e pode ser muito inclinada: o CytoNorm levou
  valores a −0,89 e a 1,55 (e a 4,96 em outra semente) num eixo de 0 a 1; o nosso nunca corrigiu
  um valor em mais de 0,17. São poucas dezenas de valores em 1,4 milhão, mas são os que aparecem
  num gráfico como eventos fora do eixo.
- **Com clusters, os dois programas diferem entre si tanto quanto cada um consigo mesmo.** O
  mapa do FlowSOM depende da semente, e um evento na fronteira de dois metaclusters recebe uma
  curva ou a outra. De uma semente para outra do mesmo programa a normalização muda de 17% a 31%
  da correção média; do nosso para o CytoNorm, de 14% a 39%. No resultado que interessa não se
  distinguem: as referências ficam a 0,013–0,014 nos dois (eram 0,022), e no simulado o erro
  contra a verdade cai de 0,085 para 0,012 (nosso) e 0,009 a 0,022 (CytoNorm, conforme o mapa).
- **Diferenças de desenho que não mudaram o resultado aqui:** o nosso treina o mapa com o mesmo
  número de eventos de cada lote (o CytoNorm reparte `nCells` igualmente entre os arquivos, mas
  usa tudo o que um arquivo menor tiver); o nosso mapa padrão é 10 × 10 e o dele 15 × 15; o
  nosso deixa sem corrigir um par (lote, cluster) com menos de 50 eventos de referência e o
  CytoNorm com 50 ou menos. No nosso
  ficaram sem corrigir de 2 a 4 pares nos dados reais e de 12 a 18 no simulado (metaclusters
  minúsculos, entre as populações), com 0,05% a 0,26% dos eventos de referência.

Nada aqui indica erro em `normalize.py`; nada foi alterado no código, só a documentação do
módulo, que agora cita estes números.

---

## 5. Autocompensação por controles de uma cor (`milc/autocomp.py`), o AutoSpill e o flowStats

Medido em 23 de novembro de 2026, com AutoSpill 0.2.0 (`carlosproca/autospill`, instalado do
GitHub com as dependências do CRAN), flowStats 4.24.0, flowCore 2.24.0 e MASS 7.3.65. Script:
`tests/reference/autocomp_compare.py` com `autocomp_autospill.R`.

### O que foi comparado

- **AutoSpill**: as duas funções do pacote, `get.marker.spillover` (regressão robusta `MASS::rlm`
  de cada detector contra o primário, nos eventos sem compensar) e `refine.spillover` (compensa,
  mede o que sobrou, corrige, repete; primeiro na escala linear até nenhuma inclinação passar de
  1e-2, depois na biexponencial até 1e-4), com os parâmetros de `get.autospill.param("minimal")`.
  Rodadas duas vezes: em todos os eventos (os algoritmos sozinhos, os mesmos eventos dos dois
  lados) e no gate de scatter do próprio AutoSpill (por tesselação), cujas linhas o R escreve e o
  Python lê, para rodar o nosso nos mesmos eventos.
- **posnegpop** do AutoSpill (`get.marker.spillover.posnegpop`): mediana dos positivos menos a
  dos negativos do próprio controle — o que o nosso `classic` faz.
- **`flowStats::spillover`** (`method = "median"`, o não marcado como base): mediana do arquivo
  marcado inteiro menos a do não marcado. É o `spillover()` que estava no flowCore e passou para
  o flowStats; o flowCore 2.24 só lê a matriz do arquivo.

### Dados

- **sim6**: seis detectores, matriz conhecida, 20.000 eventos por controle, 60% positivos,
  autofluorescência, ruído de fótons e eletrônico, 15% de debris e 1% de eventos brilhantes
  soltos, e um controle não marcado (`autocomp.simulate`).
- **sim12**: doze detectores, matriz mais densa (até 35% entre vizinhos), 15.000 eventos, limpo.
- **compdata**: os controles reais que o flowCore traz (`extdata/compdata`, um FACSCalibur,
  quatro detectores FL1-H a FL4-H, amplificação logarítmica de 4 décadas, 7.500 a 10.000
  eventos, e um tubo não marcado). Lidos pelo MIL-C, que aplica o `$PnE`, e reescritos como
  valores lineares, para que os dois programas leiam os mesmos números (o AutoSpill lê com
  `transformation = NULL` e não linearizaria).

Não houve como usar os controles do artigo do AutoSpill (FlowRepository FR-FCM-Z2SS, Z2ST, Z2SV):
a página lista os arquivos, mas baixar exige conta, e a API do FlowRepository recusa clientes não
cadastrados. Nenhuma conta foi criada.

### Concordância medida

Maior diferença entre duas matrizes, fora da diagonal, em fração (0,001 = 0,1 ponto percentual):

| | sim6 | sim12 | compdata |
|---|---|---|---|
| MIL-C AutoSpill × AutoSpill refinado, todos os eventos | 1,9e-5 | 5,5e-6 | 2,7e-4 |
| MIL-C AutoSpill × AutoSpill refinado, gate do AutoSpill | 3,0e-6 | 5,0e-6 | 2,1e-4 |
| MIL-C clássico × posnegpop do AutoSpill | 8,7e-4 | 1,3e-3 | 0,035 |
| MIL-C clássico (não marcado como negativo) × flowStats | 0,058 | 2,9e-3 | 0,017 |
| primeiro passo × refinado, ambos do AutoSpill | 6,7e-4 | 3,9e-4 | 3,2e-3 |

O número de iterações é o mesmo dos dois lados em todas as seis rodadas (5, 3, 3, 4, 5, 5). Nos
dados do FACSCalibur a diferença maior (2,7e-4) vem, até onde se pode ver, dos valores
quantizados da amplificação logarítmica: muitos eventos empatados mudam a mediana dos resíduos
que a regressão de Huber usa como escala.

Distância à matriz verdadeira (só nas simuladas):

| | sim6 | sim12 |
|---|---|---|
| MIL-C AutoSpill, todos os eventos | 0,0026 | 0,00037 |
| MIL-C AutoSpill, gate de scatter do MIL-C | 0,00017 | 0,00037 |
| AutoSpill refinado, todos os eventos | 0,0026 | 0,00037 |
| AutoSpill refinado, gate do AutoSpill | 0,00024 | 0,00036 |
| MIL-C clássico, todos os eventos | 0,0034 | 0,0013 |
| MIL-C clássico, gate do MIL-C | 0,0012 | 0,0012 |
| posnegpop do AutoSpill | 0,0030 | 0,0012 |
| flowStats | 0,058 | 0,0020 |

Em `tests/test_autocomp.py`, com outra simulação de seis detectores e sementes 0 a 7, o AutoSpill
ficou a no máximo 0,0004 da verdade e o clássico a 0,0026; o teste pede 0,001 e 0,005. O
refinamento reduz a maior inclinação que sobra nos controles compensados de 3e-4 (primeiro
passo) para 3e-5 a 1,7e-4.

Tempo, na mesma máquina e com ela ocupada: o nosso levou de 0,1 a 4 s por conjunto; o AutoSpill,
de 1 a 15 s, mais 3 a 16 s do gate por tesselação.

### Onde diferem e por quê

- **O flowStats em sim6 (5,8 pontos de erro)**: ele toma a mediana do arquivo inteiro. Com 40% de
  células negativas e 15% de debris no tubo, a mediana não é a dos positivos. É o desenho do
  método — feito para tubos quase todos positivos, como esferas —, não um defeito; no sim12, sem
  debris, fica a 0,2 ponto.
- **Clássico × posnegpop**: os dois tomam positivos e negativos do próprio controle, mas acham o
  vale de forma um pouco diferente (o nosso: duas médias e o mínimo de um histograma suavizado;
  o deles: `kmeans` e o mínimo de `density`). Nas simuladas a diferença é de 0,1 ponto. No
  compdata chega a 3,5 pontos na linha do FL3, cujo controle é fraco (mediana 63 contra 4 do não
  marcado): positivos e negativos se sobrepõem e a posição do corte pesa.
- **Clássico × AutoSpill no compdata**: o FL3 no FL4 sai 0,12 pelo AutoSpill e 0,15 a 0,23 pelo
  clássico, conforme o negativo (do próprio tubo ou o não marcado). Sem matriz verdadeira não há
  como dizer qual está certo; o critério do AutoSpill (inclinação que sobra) favorece o próprio
  AutoSpill por construção, e por isso não foi usado como juiz aqui.
- **O gate de scatter**: com 15% de debris (sim6) o gate muda tudo — o AutoSpill vai de 0,26 a
  0,02 ponto de erro. Com o nosso gate (a elipse robusta de `autogate.gaussian`, 95%) ou com o
  do AutoSpill (tesselação), as duas implementações ficam a menos de 0,03 ponto da verdade.

### O que não foi conferido

- Controles de citômetro espectral ou com mais de 12 detectores; controles em esferas.
- A fase biexponencial usa a biexponencial padrão do AutoSpill (a mesma do MIL-C) para todos os
  detectores; o AutoSpill aceita parâmetros por detector num CSV, que aqui não se lê.
- O controle de autofluorescência do AutoSpill (antígeno "AutoF") e o gate por controle em CSV
  não foram reproduzidos.
- Nenhuma matriz calculada aqui foi posta ao lado da do FlowJo ou do FACSDiva para os mesmos
  controles.

## 6. Ciclo celular (`milc/cellcycle.py`) e cinética (`milc/kinetics.py`)

### O que foi comparado

Não há aqui nenhum programa que ajuste histogramas de DNA como o ModFit ou o FlowJo, e não achei
pacote de R ou Python que faça Watson ou Dean-Jett-Fox. A referência do ciclo foi então o
próprio dado: arquivos com pulso de BrdU e iodeto de propídio, em que a fase de cada célula se lê
no gráfico BrdU × PI sem modelo nenhum. O MIL-C vê só o PI. Script:
`tests/reference/cellcycle_compare.py`.

A leitura do BrdU, escrita aqui (os autores não publicaram a deles por arquivo): singletos por
PI-H/PI-A acima de 0,82 da mediana no pico G1; BrdU+ = log10 FITC-A mais de 0,4 acima da faixa
BrdU−, traçada entre a moda do pico G1 e a do G2; BrdU+ é S, e as BrdU− abaixo ou acima do meio
dos dois picos são G1 ou G2. Isto não é a verdade do conteúdo de DNA: células que incorporaram
BrdU e chegaram a G2 (ou se dividiram) durante o pulso contam como S, e o começo de S com pouco
BrdU conta como G1. Onde um fármaco para a replicação (XL413, AZD6738, olaparibe), há células em
S sem BrdU e a referência quebra — por isso o resumo também é dado só para os controles
(siControl, DMSO, shControl).

### Dados

Yang & Hill, *The ORFIUS complex regulates ORC2 localization at replication origins*: 223
arquivos FCS de linhagens de câncer de ovário (OV8, FT194, UWB1.289 com e sem BRCA1, CaOV3,
SKOV3), BD Fortessa, 10 a 20 mil eventos cada. Dryad doi:10.5061/dryad.sn02v6x89 (CC0), também
no Zenodo (registro 10499157). O PI-A satura em 262.143 em alguns por cento dos eventos (4 e 7% nos dois arquivos olhados de perto; agregados, na maioria): esses
eventos ficam fora de todos os bins. Baixados para uma pasta própria e lidos com `python3 -I`.

### Concordância medida

Modelo menos BrdU, em pontos (média ± DP entre arquivos):

| | G1 | S | G2 | χ² reduzido (mediana) |
|---|---|---|---|---|
| **controles (39)**; BrdU: 39,8 / 39,8 / 20,4 | | | | |
| singletos, Watson | −2,8 ± 14,8 | +3,9 ± 13,5 | −1,1 ± 9,0 | 2,60 |
| singletos, Dean-Jett-Fox | +4,9 ± 4,1 | −6,6 ± 5,8 | +1,7 ± 5,2 | 1,87 |
| singletos, DJF + debris | +4,9 ± 4,0 | −6,7 ± 5,8 | +1,8 ± 5,3 | 1,41 |
| todos os eventos, DJF | +0,8 ± 4,0 | −11,4 ± 7,2 | +10,7 ± 7,6 | 9,98 |
| todos os eventos, DJF + agregados | +4,3 ± 4,2 | −6,5 ± 6,1 | +2,1 ± 5,9 | 1,63 |
| todos os eventos, DJF + debris + agregados | +4,4 ± 4,1 | −6,7 ± 6,0 | +2,3 ± 5,6 | 1,44 |
| **todos os arquivos (223)**; BrdU: 45,9 / 33,0 / 21,1 | | | | |
| singletos, Watson | −6,2 ± 16,3 | +10,9 ± 17,2 | −4,6 ± 8,8 | 2,89 |
| singletos, Dean-Jett-Fox | +2,5 ± 6,8 | −1,5 ± 11,8 | −1,1 ± 7,8 | 2,11 |
| todos os eventos, DJF + debris + agregados | +1,7 ± 6,8 | −2,6 ± 11,7 | +0,8 ± 7,5 | 1,62 |

Agregados: nos eventos sem gate, o modelo acha 24,8% (DP 5,2) de agregados na faixa ajustada; o
gate de singletos tira 25,6% (DP 4,7) da mesma faixa; r = 0,90 entre arquivos. Ajustar todos os
eventos com agregados dá, arquivo a arquivo, as fases dos singletos (S −1,1 ± 4,9; G2 +1,9 ±
5,3); sem agregados, o G2 sobe 10,4 ± 7,1 pontos — os dobletes G1+G1 contados como G2.

Cinética (`tests/reference/kinetics_compare.py`): em 23 desses arquivos (um em cada dez), a
mediana do FITC-A e a da razão FITC-A/PI-A em 100 fatias de tempo, média móvel de 5, linha de
base, pico, tempo do pico e AUC, comparados com as mesmas definições escritas de novo em R, com
os arquivos lidos pelo flowCore: maior diferença relativa 3,5e-14.

### Onde diferem e por quê

- **S abaixo do BrdU nos controles (−6,6).** É o sentido esperado da diferença entre as duas
  medidas: o BrdU marca quem passou por S durante o pulso, inclusive quem já está em G2 ou de
  volta em G1. Sem o ModFit não sei quanto dessa diferença ele também teria.
- **Watson.** Na média dos controles erra pouco, mas espalha três vezes mais. O G2 destas
  linhagens é assimétrico (cauda à esquerda); o Watson ajusta só o lado direito, faz um G2
  estreito demais e o resto vira S. O χ² dele, só dos lados ajustados, não acusa.
- **Fármacos.** Em UWB com XL413/AZD6738/olaparibe o Dean-Jett-Fox dá 17 a 18 pontos de S a
  mais que o BrdU (o Watson, 29): células paradas em S não incorporam BrdU. Aqui a referência é
  que falha.
- **Debris** quase não muda as fases destes arquivos (há pouco), mas baixa o χ² de 1,87 para
  1,41: o que sobrava abaixo do G1 era debris.

### O que não foi feito

- Nada foi posto ao lado do ModFit nem do FlowJo, que não estão disponíveis aqui.
- A população aneuploide só foi conferida em simulação: estes arquivos têm uma população cada.
- A cinética não foi conferida num fluxo de cálcio real (o FlowRepository não respondeu; o
  Zenodo não tinha um), e o R refaz as mesmas definições — escritas por mim, não por outro.

---

## O que não foi feito

- O Rphenograph não foi instalado: o pacote Python dos autores funcionou e era o preferido.
- O flowAI recebeu só as checagens de vazão e de sinal; a terceira (faixa dinâmica) não tem
  correspondente no nosso código e ficou fora.
- A normalização não foi comparada com arquivos FCS brutos e `transformList` do CytoNorm: os
  eventos foram entregues já transformados, para comparar a normalização e não a transformação.
- Os exemplos de controles do artigo do AutoSpill não foram obtidos (FlowRepository exige conta).
- Os tempos de execução foram anotados com a máquina ocupada por outros processos e não servem
  como medida de desempenho.
