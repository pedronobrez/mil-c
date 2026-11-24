---
title: Cinética
---

A cinética acompanha um parâmetro ao longo do tempo da aquisição, como num fluxo de cálcio: as células correm por um tempo, um estímulo é adicionado, e o sinal — ou a razão entre dois sinais, como no Indo-1 — sobe e desce. A aquisição é cortada em fatias de tempo, um número é tirado dos eventos de cada fatia, e a série é medida entre uma linha de base e uma resposta.

## Onde fica

Na [[gates-tab]], escolha a população, selecione as amostras na lista e depois **Model…** ▸ **Kinetics…**. Os arquivos precisam ter um parâmetro Time. O **Parameter** começa no parâmetro de X, a menos que seja o próprio Time. Os valores são os que os arquivos guardam, compensados, sem escala de exibição.

## Os ajustes

| Controle | Padrão | O que faz |
|---|---|---|
| **Parameter** | o de X | O que se acompanha. |
| **Divided by** | — | Um segundo parâmetro: acompanha-se a razão entre os dois, evento a evento (violeta sobre azul do Indo-1, por exemplo). Eventos cujo denominador não é maior que zero ficam de fora. |
| **Statistic** | Median | **Median**, **Mean** ou **% above threshold** de cada fatia. |
| **Threshold** | 0 | Para **% above threshold**: os eventos acima disto, nas unidades do arquivo (ou da razão), são os que responderam. |
| **Slice** | 100 slices | A duração de uma fatia, nas unidades do eixo. No mínimo, a aquisição é cortada em 100 fatias. |
| **Smoothing** | Moving average, 5 | **None**, **Moving average** ou **Gaussian**; o número é em fatias: quantas a média móvel abrange (arredondado para ímpar, centrado), ou o DP da gaussiana. |
| **Baseline** | o primeiro sexto | De e até: as fatias antes do estímulo. |
| **Response** | o resto | De e até: as fatias depois dele. |

As faixas de **Baseline** e **Response** também podem ser arrastadas no gráfico. Todas as amostras selecionadas são cortadas nas mesmas fatias, no intervalo que cobre todas.

O tempo fica em segundos quando todos os arquivos dizem quanto dura a sua unidade de Time (`$TIMESTEP`); senão, nas unidades dos próprios arquivos, e o eixo diz isso.

## A série

Em cada fatia, a estatística é tirada dos eventos da população. Uma fatia com menos de 5 eventos fica sem número. A série é então suavizada, cada fatia com o mesmo peso e as fatias vazias puladas. Os pontos do gráfico são as fatias, a linha é a série suavizada, e cada amostra tem a sua cor (ou cada valor da coluna **Colour by** da aba Gates).

Cada fatia tem um erro-padrão: DP/√n para a média, 1,2533 × DP robusto/√n para a mediana (exato para a distribuição normal, uma indicação nas outras), √(p(1−p)/n) para uma porcentagem. O erro de uma fatia suavizada sai dos erros das fatias que entram nela.

## A tabela

| Coluna | Significado |
|---|---|
| **Events** | Eventos da população acompanhada. |
| **Baseline**, **Baseline SD** | Média e DP das fatias da linha de base, sem suavizar: quanto uma fatia varia sozinha. |
| **Peak**, **Peak SE** | O maior valor da série suavizada na resposta, com o seu erro-padrão. |
| **Peak − baseline**, **Peak / baseline** | A subida. |
| **Time of peak** | No eixo do gráfico. |
| **Time to peak** | Contado do início da resposta, onde está o estímulo. |
| **Slope**, **Slope SE** | De uma reta pelas fatias (sem suavizar) do início da resposta até o pico, por mínimos quadrados, por unidade de tempo. Precisa de três fatias. |
| **AUC** | Área entre a série suavizada e a linha de base ao longo da resposta (trapézios): sinal × tempo. |

As suas colunas da [[samples-tab]] vêm primeiro, para que a tabela possa ir direto para **Compare groups…** (que começa em **Peak − baseline**; veja [[statistics-tab]]), **Export CSV…** ou **Copy** ([[tables-and-prism]]). **Figure…** faz do gráfico uma figura, com as faixas ([[figures]]).

## Como foi conferido

Em 23 arquivos públicos com parâmetro Time, a mediana de um parâmetro e a de uma razão entre dois em 100 fatias, média móvel de 5, linha de base, pico, tempo do pico e AUC foram calculados de novo em R a partir das mesmas definições, com os arquivos lidos pelo flowCore; a maior diferença relativa foi 3,5 × 10⁻¹⁴ (`docs/comparacao.md`, seção 6). Isso confere a aritmética, não um fluxo de cálcio de verdade: nenhum dado público de fluxo de cálcio pôde ser obtido, e o código em R repete definições escritas aqui.
