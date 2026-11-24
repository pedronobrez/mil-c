---
title: Proliferação
---

A proliferação é lida da diluição de um corante (CFSE, CellTrace Violet e semelhantes). Uma célula que se divide dá a cada filha metade do seu corante, então as gerações são picos a um fator de dois um do outro: igualmente espaçados num eixo logarítmico. Um modelo do histograma diz quantas células há em cada geração, e delas saem os índices do FlowJo.

## Onde fica

Na [[gates-tab]], escolha a população (em geral, células únicas vivas), selecione as amostras na lista e depois **Model…** ▸ **Proliferation…**. O diálogo ajusta todas as amostras selecionadas e mostra uma por vez. O **Parameter** começa no parâmetro de X. Os valores são os que os arquivos guardam, compensados, sem escala de exibição; o modelo trabalha no logaritmo deles, então eventos em zero ou abaixo ficam de fora.

## Os ajustes

| Controle | Padrão | O que faz |
|---|---|---|
| **Parameter** | o de X | O corante. |
| **Undivided peak** | Fitted in each sample | Uma amostra que não foi estimulada: o pico dela diz qual pico das outras é o não dividido. |
| **Generations, at most** | 8 | Picos no modelo depois do não dividido, de 1 a 12. |
| **Autofluorescence** | 0 | O sinal de células sem o corante, nas unidades do arquivo: só o corante acima dele se divide pela metade. 0: nenhum. |
| **Peaks widen with each generation** | desligado | Desmarcado, todo pico tem a mesma largura no eixo logarítmico. Marcado, a largura pode crescer na mesma fração, até 30%, a cada geração. |

## O modelo

O histograma do log10 do corante (256 classes) é ajustado como uma soma de gaussianas, uma para cada geração *g* = 0, 1, 2…, em

`posição(g) = posição(0) − g · log10(2)`

todas com a mesma largura (ou com uma que cresce com *g*, quando marcado). Com **Autofluorescence** acima de zero, o que se divide pela metade é o corante acima dela: `posição(g) = log10(fundo + (10^posição(0) − fundo) / 2^g)`. O número de eventos de cada geração é achado por mínimos quadrados não negativos; só a posição da geração 0, a largura e o seu crescimento são procurados.

Saber qual pico é a geração 0 é a parte difícil:

- **Com um controle** escolhido em **Undivided peak**, uma gaussiana é ajustada ao histograma do controle, e a geração 0 de cada amostra é procurada a até um terço de divisão dela (log10(2) / 3). Células estimuladas que não se dividiram raramente têm exatamente o brilho das que ficaram sem estímulo, mas nunca estão uma divisão inteira longe. A largura é sempre a da própria amostra. O pico do controle aparece no gráfico como uma linha tracejada.
- **Sem controle**, toda posição é tentada e, das que ajustam tão bem quanto a melhor (a 5% dela), fica a mais fraca: o pico mais brilhante que contém células é tomado como o não dividido. Isso erra quando todas as células se dividiram. Use um controle sempre que houver.

## A tabela

Uma linha por amostra, com as suas colunas da [[samples-tab]]. *N(g)* é o número de eventos na geração *g*; as células que existiam no começo e das quais esses descendem são *N(g) / 2^g*.

| Coluna | Significado | Fórmula |
|---|---|---|
| **Generations** | A geração mais alta com pelo menos 0,5% dos eventos modelados | |
| **% Divided** | Fração das células do começo que se dividiram pelo menos uma vez (frequência de precursores) | 100 · Σ(g≥1) N/2^g ÷ Σ(g≥0) N/2^g |
| **Division Index** | Divisões por célula do começo, incluídas as que nunca se dividiram | Σ(g≥1) g·N/2^g ÷ Σ(g≥0) N/2^g |
| **Proliferation Index** | Divisões por célula que se dividiu | Σ(g≥1) g·N/2^g ÷ Σ(g≥1) N/2^g |
| **Expansion Index** | Quantas vezes a cultura cresceu | Σ(g≥0) N ÷ Σ(g≥0) N/2^g |
| **Replication Index** | Quantas vezes as células que se dividiram cresceram | Σ(g≥1) N ÷ Σ(g≥1) N/2^g |
| **Undivided peak** | Onde está a geração 0, nas unidades do arquivo | |
| **Peak SD (decades)** | A largura de um pico | |
| **RMS**, **RMS %** | O resíduo, em eventos por classe e como % da classe mais alta | |
| **% modelled** | Fração dos eventos que as gerações explicam | |
| **% in generation 0**, **1**… | Fração dos eventos modelados em cada geração | |

São as definições do FlowJo. O flowFit do R chama o Expansion Index de "proliferation index". Cada índice vem com um erro-padrão (**SE**) de contagem: quanto as mesmas células, adquiridas de novo, mudariam o número. Ele não diz se os picos são gaussianos; o RMS diz. Células que morreram não são vistas, então todo índice é das células que estão lá para serem contadas.

## O gráfico e os botões

O gráfico mostra o histograma num eixo logarítmico, cada geração e o modelo, para a amostra escolhida na tabela; a linha acima dá as gerações, a % dividida, os dois índices principais, o RMS e onde está o pico não dividido. **Figure…** faz dele uma figura ([[figures]]); **Compare groups…** testa uma coluna (de início o **Division Index**) entre grupos de amostras ([[statistics-tab]]); **Export CSV…** e **Copy** levam a tabela para fora ([[tables-and-prism]]).

## Limites

O modelo supõe picos gaussianos de uma só largura no eixo log, o que vale para um corante bem marcado e compensado. Um corante que vaza para um canal vizinho, ou uma população com o pico não dividido largo, dá picos que não são gaussianos: olhe o RMS e o gráfico. Sem controle, uma cultura em que todas as células se dividiram é lida com uma geração de diferença. Nenhuma comparação com a plataforma de proliferação do próprio FlowJo foi feita.
