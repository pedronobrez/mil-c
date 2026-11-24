---
title: Glossário
---

Os termos que você encontra no MIL-C, na ordem alfabética do termo em inglês (entre parênteses quando o termo é dado em português), cada um com a página onde é usado.

## A – C

- **Adjusted p (p ajustado)** — um valor de p corrigido pelo número de testes feitos ao mesmo tempo, para que testar vinte clusters não produza um resultado "significativo" só por acaso. Ver [[cluster-statistics]].
- **ANOVA** — análise de variância: um teste de se as médias de três ou mais grupos diferem. A ANOVA de Welch admite variâncias diferentes; a de medidas repetidas é para os mesmos sujeitos em todos os grupos. Ver [[cluster-statistics]].
- **Arcsinh** — o seno hiperbólico inverso, asinh(x / cofator): linear perto de zero, logarítmico para valores grandes, definido para valores negativos. Uma das escalas de exibição e de análise. Ver [[scales-and-transforms]].
- **Autocompensação (autocompensation)** — calcular a matriz de spillover a partir de controles de uma cor em vez de digitá-la. Ver [[autocompensation]].
- **Autofluorescência (autofluorescence)** — a luz que as células emitem sem corante; no unmixing espectral ela pode ser mais um espectro. Ver [[spectral-unmixing]].
- **Backgating** — mostrar onde os eventos de uma população caem em cada gate acima dela. Ver [[gates-tab]].
- **Normalização de lotes (batch normalisation)** — trazer amostras adquiridas em dias diferentes (lotes) umas sobre as outras, aprendendo com amostras de referência que deveriam ser iguais em todo lote. Ver [[normalisation]].
- **Benjamini-Hochberg** — o ajuste que controla a taxa de falsas descobertas; o padrão no MIL-C. Ver [[cluster-statistics]].
- **Biexponencial (biexponential)** — a escala de exibição do FlowJo, linear em torno de zero e logarítmica acima, com uma largura e um número de décadas. A escala padrão do MIL-C para marcadores, e a que os algoritmos recebem, como números de 0 a 1023. Ver [[scales-and-transforms]].
- **Bonferroni, Holm** — ajustes que controlam a chance de qualquer falso positivo; mais estritos que o Benjamini-Hochberg. Ver [[cluster-statistics]].
- **Gate booleano (boolean gate)** — uma população feita de outras: os eventos em todas, em alguma ou em nenhuma delas. Ver [[drawing-gates]].
- **Ciclo celular (cell cycle)** — G1, S e G2 ajustados a um histograma de conteúdo de DNA. Ver [[cell-cycle]].
- **Cluster** — um grupo de eventos parecidos nos marcadores, achado por um algoritmo (FlowSOM, PhenoGraph) ou desenhado à mão. Ver [[clustering-methods]] e [[clusters-tab]].
- **Cofactor (cofator)** — o divisor em arcsinh(x / cofator), que define onde a escala passa de linear a logarítmica: 150 é o usual para fluorescência, 5 para citometria de massa. Ver [[scales-and-transforms]].
- **Compensação (compensation)** — retirar de cada canal de fluorescência a luz que vazou para ele de outros fluorocromos, pela matriz de spillover. Os canais compensados se chamam `Comp-…`. Ver [[compensation]].
- **Concatenação (concatenation)** — juntar os eventos de muitas amostras num conjunto só, cada evento lembrando a sua amostra, para as reduções de dimensão e clusterizações. Ver [[concatenation]].
- **Clusterização de consenso (consensus clustering)** — o modo do FlowSOM de juntar os nós do seu mapa em metaclusters: muitas clusterizações hierárquicas reamostradas, combinadas. Ver [[clustering-methods]].
- **Count (contagem)** — o número de eventos numa população. Ver [[statistics-tab]].
- **CV** — coeficiente de variação: 100 × desvio-padrão / média. O CV robusto é 100 × desvio-padrão robusto / mediana. Ver [[statistics-tab]].

## D – H

- **Abundância diferencial (differential abundance)** — se um cluster é uma fração maior ou menor da amostra num grupo do que noutro, testado sobre contagens. Ver [[abundance]].
- **Expressão diferencial (differential expression)** — se as células de um cluster expressam mais ou menos de um marcador num grupo do que noutro, testado sobre a mediana de cada amostra. Ver [[differential-expression]].
- **Dispersão (dispersion)** — num modelo binomial negativo, quanto as contagens variam entre amostras além do que o sorteio de células explicaria. Ver [[abundance]].
- **Downsampling (subamostragem)** — tirar menos eventos do que uma amostra tem. No MIL-C isso é feito pelo plano de amostragem, que balanceia as amostras. Ver [[concatenation]].
- **Embedding (redução de dimensão)** — um mapa em duas dimensões dos eventos, como o UMAP ou o tSNE. Ver [[embeddings]].
- **Evento (event)** — uma partícula que passou pelo laser: em geral uma célula. Uma linha de um arquivo FCS. Ver [[opening-data]].
- **FCS** — Flow Cytometry Standard, o formato de arquivo que os citômetros gravam. O MIL-C lê as versões 2.0, 3.0 e 3.1; os seus arquivos são só lidos, e as exportações são arquivos novos. Ver [[opening-data]].
- **FDR** — taxa de falsas descobertas: entre os resultados chamados significativos, a fração esperada que não é. O p ajustado por Benjamini-Hochberg costuma ser chamado de FDR ou valor q. Ver [[cluster-statistics]].
- **FlowSOM** — clusterização por um mapa auto-organizável cujos nós são juntados num número escolhido de metaclusters. Ver [[clustering-methods]].
- **FMO** — fluorescence minus one: um controle corado com todos os fluorocromos menos um, para ver onde terminam os negativos daquele. O MIL-C não tem ferramenta específica para FMO; sobreponha o controle a uma amostra em Gates para posicionar um gate. Ver [[gates-tab]].
- **Fold change (razão)** — a razão entre grupos; na abundância diferencial é dada em log2, então 1 é o dobro da fração e −1 a metade. Ver [[abundance]].
- **Gate** — uma região em um ou dois parâmetros que seleciona eventos; gates uns dentro dos outros formam a árvore de gates. Ver [[drawing-gates]] e [[auto-gates]].
- **Gating-ML** — o formato padrão de arquivo para gates, lido pelo FlowKit, pelo flowCore e pelo Cytobank. Ver [[workspaces-and-gating-ml]].
- **Geometric mean (média geométrica)** — exp da média dos logaritmos; eventos em zero ou abaixo ficam de fora. Ver [[statistics-tab]].
- **Greenhouse-Geisser** — a correção de uma ANOVA de medidas repetidas para grupos cujas diferenças não se espalham igualmente. Ver [[cluster-statistics]].
- **Group (grupo)** — as amostras que compartilham um valor de uma coluna (ou de várias) na tabela de Samples. Ver [[samples-tab]] e [[sample-order]].

## J – M

- **Jaccard** — a fração que dois conjuntos têm em comum. O PhenoGraph pesa cada ligação pelo índice de Jaccard das vizinhanças dos dois eventos. Ver [[clustering-methods]].
- **Cinética (kinetics)** — um parâmetro acompanhado ao longo do tempo, como num fluxo de cálcio. Ver [[kinetics]].
- **Kruskal-Wallis, Mann-Whitney** — testes de postos, para três ou mais grupos e para dois, que não supõem nada sobre a forma dos valores. Ver [[cluster-statistics]].
- **Leiden, Louvain** — algoritmos que cortam um grafo em comunidades de alta modularidade; o último passo do PhenoGraph. O Leiden é o padrão. Ver [[clustering-methods]].
- **Logicle** — uma escala do tipo biexponencial com parâmetros padronizados (W, M, A, T). Ver [[scales-and-transforms]].
- **MAD** — desvio absoluto mediano: a mediana das distâncias à mediana. Uma estatística da tabela, e a unidade dos limites do controle de qualidade. Ver [[statistics-tab]] e [[quality-control]].
- **Median (mediana)** — o valor do meio dos eventos: o número usual para o nível de um marcador, robusto a valores extremos. Ver [[statistics-tab]].
- **Metacluster** — um cluster do FlowSOM: um grupo de nós do mapa. Na aba Clusters, simplesmente um cluster. Ver [[clustering-methods]].
- **MFI** — intensidade média (ou mediana) de fluorescência. O MIL-C relata mediana, média ou média geométrica pelo nome, em vez de "MFI". Ver [[statistics-tab]].
- **t moderado (moderated t)** — a estatística t do limma, cuja variância é puxada em direção ao que todas as linhas sugerem; vale graus de liberdade a mais com poucas amostras. Ver [[differential-expression]].
- **Modularidade (modularity)** — quanto mais peso uma partição de um grafo mantém dentro das suas comunidades do que o acaso manteria. O PhenoGraph a maximiza. Ver [[clustering-methods]].

## N – R

- **Binomial negativa (negative binomial)** — a distribuição de contagens usada na abundância diferencial: parecida com a de Poisson, com variação extra entre amostras. Ver [[abundance]].
- **Preservação de vizinhos (neighbour preservation)** — a fração dos 15 vizinhos mais próximos de cada evento que continuam vizinhos num mapa: como o MIL-C avalia uma redução de dimensão. Ver [[embeddings]].
- **Desenho pareado (paired design)** — os mesmos animais ou doadores medidos em todos os grupos; os testes olham então a mudança de cada um. Ver [[cluster-statistics]].
- **PCA** — análise de componentes principais: as direções em que os eventos mais variam, cada uma uma mistura legível de marcadores. Ver [[embeddings]].
- **Perplexity (perplexidade)** — a configuração do tSNE, aproximadamente o número de vizinhos que cada evento considera. Ver [[embeddings]].
- **PhenoGraph** — clusterização sem escolher o número: um grafo de vizinhos mais próximos, pesos de Jaccard e comunidades. Ver [[clustering-methods]].
- **Plate (prancha)** — vários painéis de figura numa página. Ver [[the-plate]].
- **Population (população)** — os eventos que um gate seleciona, com a sua contagem e porcentagens. Ver [[gates-tab]].
- **Prism, .pzfx** — o GraphPad Prism e o seu formato de arquivo, que o MIL-C grava. Ver [[tables-and-prism]].
- **Proliferação (proliferation)** — gerações de um corante diluído a cada divisão, ajustadas a um histograma. Ver [[proliferation]].
- **Controle de qualidade (quality control)** — achar trechos de uma aquisição em que a vazão ou o sinal deram errado, e deixá-los de fora. Ver [[quality-control]].
- **Resolution (resolução)** — a configuração do PhenoGraph: acima de 1, mais clusters e menores; abaixo, menos e maiores. Ver [[clustering-methods]].
- **Robust SD (desvio-padrão robusto)** — metade da distância entre os percentis 15,87 e 84,13: o desvio-padrão de uma normal com essa largura, surdo a valores extremos. Ver [[statistics-tab]].

## S – W

- **Plano de amostragem (sampling plan)** — quantos eventos cada amostra cede à concatenação, e por quê. Ver [[concatenation]].
- **Seed (semente)** — o número que fixa as escolhas aleatórias de um algoritmo, de modo que a mesma semente dá o mesmo resultado. Ver [[explore-tab]].
- **Mapa auto-organizável (self-organising map)** — uma grade de nós treinada para cobrir os eventos, o primeiro passo do FlowSOM. Ver [[clustering-methods]].
- **Controle de uma cor (single-stain control)** — uma amostra corada com um único fluorocromo, de onde se mede o spillover ou um espectro. Ver [[autocompensation]] e [[spectral-unmixing]].
- **Spillover** — a fração da luz de um fluorocromo lida no detector de outro. A matriz de spillover a guarda para cada par. Ver [[compensation]].
- **Unmixing espectral (spectral unmixing)** — num citômetro espectral, recuperar o sinal de cada fluorocromo a partir de muitos detectores e dos seus espectros de referência. Ver [[spectral-unmixing]].
- **SD** — desvio-padrão, com n − 1. Ver [[statistics-tab]].
- **tSNE** — t-distributed stochastic neighbour embedding, um mapa que preserva as vizinhanças locais. Ver [[embeddings]].
- **UMAP** — uniform manifold approximation and projection, um mapa montado a partir de um grafo de vizinhos mais próximos. Ver [[embeddings]].
- **Volcano** — razão ou diferença na horizontal, −log10 do p ajustado na vertical; os clusters que importam ficam nos cantos de cima. Ver [[abundance]] e [[differential-expression]].
- **Teste t de Welch (Welch's t test)** — um teste t que não supõe que os dois grupos tenham variâncias iguais; o padrão do MIL-C para dois grupos. Ver [[cluster-statistics]].
- **Workspace** — um arquivo `.wsp` do FlowJo: amostras, compensação, escalas e gates. O MIL-C os lê e os grava. Ver [[workspaces-and-gating-ml]].
