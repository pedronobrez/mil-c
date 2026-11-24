---
title: A aba Explore
---

A aba Explore roda os algoritmos sobre os eventos concatenados e mostra o resultado: quaisquer dois eixos, coloridos por grupo, amostra, cluster ou expressão de um marcador, num painel só ou num por grupo. Enquanto as amostras não forem concatenadas na [[samples-tab|aba Samples]] (ver [[concatenation]]) ela avisa, e os seus botões ficam desligados.

A aba tem três partes: as configurações à esquerda, o gráfico no meio e o modo de olhar à direita.

## Esquerda: o que calcular

### Markers

A lista dos canais de fluorescência corados — dispersão e tempo não são marcadores. Os que têm nome de marcador vêm marcados. São esses os marcadores que os algoritmos veem: desmarque o que foi usado para fazer gates (viabilidade, linhagem), que de outro modo dominaria o mapa. Pelo menos dois precisam estar marcados.

### Escala e semente

| Controle | Opções | Significado |
|---|---|---|
| **Scale** | **Biexponential** (padrão), **Arcsinh, cofactor 150**, **Logicle, W 0.5, M 4.5** | Como os valores são transformados antes de os algoritmos os verem. Biexponential é cada marcador na escala que o seu eixo tem em Gates — a biexponencial do FlowJo, a menos que outra tenha sido escolhida — como números de 0 a 1023, que é o que o FlowJo entrega aos seus plugins. Logicle é a Logicle padrão em 0 a 1023. Ver [[scales-and-transforms]] |
| **Seed** | 0 a 999.999, padrão 0 | A semente aleatória. A mesma semente, marcadores e configurações dão o mesmo resultado |

Nada mais é reescalado: os marcadores já chegam numa escala única.

### Os algoritmos

Cada um tem as suas configurações e o seu botão **Run**. Eles rodam em segundo plano, e vários podem rodar ao mesmo tempo; enquanto isso o botão diz "Running…".

| Algoritmo | Configurações (padrão) | Dá |
|---|---|---|
| UMAP | **Neighbours** (15), **Min. distance** (0.5) | Dois eixos, UMAP1 e UMAP2 |
| tSNE | **Perplexity** (30) | Dois eixos, tSNE1 e tSNE2 |
| FlowSOM | **Clusters** (5), **Grid** (10: um mapa de 10 × 10) | Um cluster para cada evento |
| PhenoGraph | **Neighbours** (30), **Resolution** (1.0), **Method** (Leiden ou Louvain), **Min. cluster size** (10) | Um cluster para cada evento, sem que se diga quantos |
| PCA | nenhuma | Dois eixos, PCA1 e PCA2, cada um com a sua fração da variância |

**Run all three** roda UMAP, tSNE e FlowSOM juntos. O que cada algoritmo faz, e o que as suas configurações mudam, está em [[embeddings]] e [[clustering-methods]]. Rodar um de novo substitui o resultado anterior; a semente, os marcadores, a escala e as configurações ficam guardados com o resultado e são salvos com o arquivo concatenado.

Abaixo do botão, algumas linhas descrevem o que foi calculado: o número de eventos e de amostras; para UMAP e tSNE, a fração dos 15 vizinhos mais próximos de cada evento que o mapa mantém; para o PhenoGraph, o número de clusters e a modularidade; para o PCA, a fração da variância de cada componente e os três marcadores que mais pesam nele, por exemplo `PCA1, 41%: +0.72 CD4, −0.55 CD8, +0.30 CD25`.

## Meio: o gráfico

Cada evento concatenado é um ponto. Uma nova redução de dimensão é mostrada assim que fica pronta; uma nova clusterização colore o gráfico.

## Direita: como olhar

| Controle | O que faz |
|---|---|
| **X**, **Y** | Os eixos: UMAP1/2, tSNE1/2, PCA1/2 (com a fração da variância no rótulo) ou qualquer parâmetro. Um marcador é desenhado na sua escala de Gates, com marcas em décadas |
| **Colour by** | **One colour**, o agrupamento escolhido em Samples, qualquer coluna das amostras que varie entre elas, **Sample**, **FlowSOM cluster**, **PhenoGraph cluster** ou — depois de um separador — a expressão de qualquer marcador |
| **Split by** | **One panel**, ou um painel por valor de qualquer das mesmas colunas, amostras ou clusters. Em cada painel os outros eventos aparecem esmaecidos atrás, para que se leia onde um grupo está contra onde todos estão |
| **and by** | Uma segunda subdivisão: um painel para cada combinação, `Treated · 7`. As linhas são os valores da primeira, as colunas os da segunda |
| **Ramp** | A rampa de cores da expressão |
| **Point size** | 1 a 14, padrão 4 |
| **Opacity** | 5% a 100%, padrão 80% |

A **Legend** lista cada cor com o seu número de eventos; para a expressão, dá a faixa desenhada. As cores da expressão vão do percentil 1 ao 99, para que poucos eventos muito brilhantes não apaguem o resto.

**Figure…** transforma o gráfico, como está, numa figura — tamanho, fonte, cores, legenda — em SVG, PDF, PNG ou TIFF: ver [[figures]].

A ordem em que os grupos aparecem na legenda e nos painéis é definida em [[sample-order]].

## Para onde vão os resultados

Os clusters aparecem na [[clusters-tab|aba Clusters]], com o seu heatmap, frequências e testes. As reduções de dimensão são oferecidas lá como mapa. A [[command-line|linha de comando]] faz o mesmo com `milc analyze`.
