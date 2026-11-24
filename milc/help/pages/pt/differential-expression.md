---
title: Expressão diferencial
---

A expressão diferencial pergunta se as células de um cluster carregam mais ou menos de um marcador num grupo do que noutro. É a companheira da [[abundance|abundância diferencial]], que pergunta se há mais ou menos dessas células. O método é o DS-limma do diffcyt (Weber et al., *Commun Biol* 2019): a mediana de cada marcador em cada cluster, por amostra, comparada entre grupos por um modelo linear e pelo t moderado do limma (Smyth, 2004).

Abra com **Differential expression…** na [[clusters-tab|aba Clusters]]. Ele usa os clusters como a aba os mostra e o agrupamento escolhido lá em **Compare** e **within each**.

## Os controles

A primeira linha:

| Controle | O que faz |
|---|---|
| **Compare** *coluna* | O grupo de que é a diferença. Uma diferença positiva é mais alta neste grupo |
| **against** | O grupo de referência |
| **Adjust for** | Opcional: outra coluna das amostras — um lote, um dia. O seu efeito é estimado ao lado do do grupo e retirado. A coluna que diz de qual animal ou doador é cada amostra torna a comparação pareada |
| **Markers** | Um menu de caixas de marcação, com o rótulo "3 of 8" e assim por diante. Todo marcador corado é oferecido; os que a clusterização não usou vêm marcados (todos, se ela usou todos) |
| **Scale** | **Biexponential** ou **arcsinh (cofactor 150)**: a escala em que as medianas são tomadas e em que fica a diferença. Começa na escala em que a clusterização rodou |

A segunda linha:

| Controle | Padrão | O que faz |
|---|---|---|
| **Cells in a sample, at least** | 10 | Uma amostra com menos células do cluster do que isto não dá mediana para ele: a mediana de três células é ruído |
| **Samples in a group, at least** | 2 | Um cluster com menos amostras assim do que isto em algum dos grupos não é testado, e a nota diz isso |
| **Weigh by cells** | desligado | Cada mediana pesa como o número de células de que ela é, como no diffcyt |
| **Correction** | Benjamini-Hochberg | Sobre todos os clusters e marcadores testados: Benjamini-Hochberg, Holm, Bonferroni ou None |

Toda mudança recalcula na hora. Com uma coluna em **within each** na aba Clusters, a análise — priori e tudo — é feita separadamente dentro de cada um dos seus valores, e os p são ajustados sobre os testes de todos os painéis juntos.

## O método

1. Para cada amostra, cluster e marcador: a mediana do marcador entre as células da amostra no cluster, na escala escolhida.
2. Para cada cluster e marcador, um modelo linear `mediana ~ intercepto + grupo (+ o fator ajustado)` sobre as amostras dos dois grupos. O coeficiente do grupo é a **diferença** de mediana entre os grupos, nas unidades da escala.
3. A variância residual de cada linha é puxada em direção ao que todas as linhas sugerem (o Bayes empírico do limma: `lmFit` + `eBayes`, sem tendência, sem a opção robusta). Isso vale graus de liberdade a mais quando há poucas amostras e muitas linhas. O t moderado dá então o p.
4. Os valores de p são ajustados sobre todos os clusters e marcadores testados juntos.

O fator ajustado é um efeito fixo: uma coluna para cada um dos seus níveis. Com o animal como fator e todo animal nos dois grupos, o teste sem moderação é o teste t pareado.

### Como ler a diferença

- **Biexponential**: as medianas estão na escala dos gráficos, como números de 0 a 1023; onde a escala é logarítmica, 100 é cerca de meia década.
- **arcsinh**: asinh(x / 150); onde o marcador é brilhante, 0,69 é uma duplicação.

## O resultado

A tabela da esquerda, uma linha por cluster e marcador:

| Coluna | Significado |
|---|---|
| *within each* | O painel, quando há um |
| Cluster, Marker | O que foi testado |
| *Grupo A* (median), *Grupo B* (median) | A média das medianas das amostras de cada grupo |
| Difference | O coeficiente do grupo: o primeiro menos o segundo, ajustado pelo resto |
| p | O p do t moderado |
| p adjusted (*método*) | Depois da correção escolhida |
| (sem nome) | Estrelas pelo p ajustado |
| n | Amostras com mediana em cada grupo, `5 vs 4` |
| Moderated t | A estatística |
| df | Graus de liberdade residuais mais os que as outras linhas valeram |
| Note | Por que uma linha não foi testada |

À direita, duas abas:

- **Heatmap** — clusters nas linhas, marcadores nas colunas, cada célula a diferença: vermelho mais alto no primeiro grupo, azul no segundo, branco nenhuma, cinza não testado. Estrelas onde o p ajustado fica abaixo de 0,05 (`*` < 0,05, `**` < 0,01, `***` < 0,001, `****` < 0,0001). A linha embaixo diz quanto vale a cor mais forte.
- **Volcano** — cada cluster e marcador um ponto, a diferença na horizontal e −log10 do p ajustado na vertical, com uma linha-guia na diferença zero; os pontos significativos são coloridos e nomeados `cluster · marcador`.

**Figure…**, no canto das abas, faz uma figura da que está à mostra (ver [[figures]]). Um parágrafo embaixo diz o que foi feito e o que não se pode concluir. **Export CSV…** (como `differential_expression.csv`), **Copy** e **Close** ficam no rodapé.

## O que ela não pode dizer

- Uma mediana é do cluster como ele foi desenhado. Um marcador sobre o qual os clusters foram construídos fica, por construção, numa faixa estreita dentro de cada cluster: o que um grupo faz com ele aparece como células mudando de cluster (abundância), e só o resto aparece aqui. Prefira marcadores que a clusterização não usou.
- Uma minoria do cluster ligando um marcador — uma mudança bimodal — mexe pouco ou nada na mediana.
- As células são as que foram concatenadas, não todas as adquiridas.

## Quão bem calibrada ela é

O t moderado coincide com o do limma 3.68 em cerca de sete algarismos significativos (ver [[measured-facts]]). Em simulações com cinco ou três amostras por grupo, sem pesos, p < 0,05 apareceu com a frequência certa quando não havia diferença (5,1% e 5,4%); com pesos, cerca de 6,5%, vezes demais, e por isso **Weigh by cells** vem desligado. Com duas amostras por grupo quase nada foi achado. O diffcyt em si não foi rodado, e o estimador mais novo do limma para linhas com amostras faltando não está implementado (aqui se segue `eBayes(legacy = TRUE)`).
