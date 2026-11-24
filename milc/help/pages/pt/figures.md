---
title: Figuras
---

Todo gráfico que importa pode virar uma figura para artigo: tamanho em milímetros, fonte, textos, cores e legenda, gravada em SVG, PDF, PNG ou TIFF. Os gráficos da janela são para olhar; a figura é desenhada de novo para impressão, a partir dos mesmos dados, então o que se exporta é o que estava na tela, com outro acabamento. Em SVG e PDF o texto continua texto, e um rótulo ainda pode ser editado no Illustrator ou no Inkscape.

## Onde ficam os botões Figure…

| Onde | O que vira figura |
|---|---|
| [[gates-tab]] | O gráfico mostrado: pseudocor, contornos, histogramas, sobreposições, violinos e caixas, com os gates desenhados. |
| **Matrix…** e **Backgate…** na aba Gates | A grade de todos os pares; a população em cada gate acima dela. |
| Diálogos **Model…** | O ajuste de [[proliferation]], [[cell-cycle]] ou [[kinetics]]. |
| [[explore-tab]] | O mapa, colorido como está na tela. |
| [[clusters-tab]] | O heatmap, as frequências por grupo, a árvore. |
| [[abundance]] e [[differential-expression]] | Os gráficos de vulcão e o heatmap de diferenças. |
| [[normalisation]] | **Figure: before…** e **Figure: after…**. |

## O diálogo Figure

A prévia à esquerda é desenhada pelo mesmo código que grava o arquivo: o que se vê é o que se grava. Embaixo, uma linha dá o tamanho e quantos pixels terá um PNG ou TIFF. Os controles que a figura não usa (uma escala de cores num histograma, por exemplo) não aparecem.

| Controle | O que faz |
|---|---|
| **Preset** | **Last used**, **Defaults**, ou um visual salvo com um nome. **Save as…** guarda o visual atual com um nome; **Delete** esquece um. Um preset leva o visual, não os textos desta figura (título, rótulos dos eixos, título da legenda). |
| **Fit** | Um tamanho: **Single column, 85 × 65 mm**, **1.5 columns, 120 × 90 mm**, **Double column, 180 × 110 mm**, **Square, 85 × 85 mm**, ou **Custom**. |
| **Width × height** | Em milímetros, de 20 a 600. |
| **Resolution** | Para PNG e TIFF, de 72 a 1200 dpi (300 por padrão). Um arquivo vetorial não tem resolução. |
| **Font**, **Size** | Qualquer fonte do sistema (Arial por padrão), em pontos (8 por padrão). |
| **Title**, **X label**, **Y label** | Vazios: os rótulos do próprio gráfico. |
| **Frame** | **Left and bottom**, **Box** ou **None**; com **Ticks** e **Grid**. |
| **Point size**, **Opacity** | Dos pontos de um mapa; a opacidade é a das barras num gráfico de frequências. |
| **Points as a bitmap inside vector files** | Para mapas com muitíssimos eventos: texto e eixos continuam vetoriais, a nuvem de pontos vira imagem, e o arquivo ainda abre. |
| **Legend** | **Right**, **Below**, **Inside** ou **Hidden**; **Counts** acrescenta o número de eventos de cada categoria. |
| **Legend title** | Vazio: o do próprio gráfico. |
| **Colour ramp** | Para o que é colorido por um valor. Começa na escala de cores que a janela mostra. |
| **Numbers in the cells** | Num heatmap. |
| **Gates, with their names and shares** | Num gráfico da aba Gates, descrito abaixo. |
| **Transparent background** | Sem branco atrás da figura. |
| Categorias | Cada categoria com a sua cor e o seu rótulo: clique duas vezes numa cor para mudá-la, digite sobre um rótulo para renomeá-lo só na figura. |

**Export…** pede um arquivo; o formato vem da extensão ou do tipo escolhido no diálogo de arquivo: **SVG**, **PDF**, **PNG** ou **TIFF**. **Close** guarda o visual como o último usado. **Add to plate**, com a lista ao lado, põe a figura como está numa prancha: uma que já existe, ou **New plate** (veja [[the-plate]]).

O visual é lembrado de uma figura para a outra: quem exporta uma figura para um artigo está prestes a exportar seis que precisam combinar. Uma figura nova começa do último visual, com os seus textos vazios e a escala de cores que a janela mostra.

## Gates numa figura

Uma figura da aba Gates leva os gates da população mostrada, desenhados como o gráfico os desenha: cada um em preto, com o nome e a fração da população (`CD4+  41.2%`) sobre um fundo claro, no mesmo lugar que na tela. Cada um dos quatro cantos de um quadrante tem o seu nome e fração. Num histograma, as faixas são desenhadas e rotuladas no alto. Os gates nunca mudam os limites dos eixos. Desmarque **Gates, with their names and shares** para ter só o gráfico.

Só os gates desenhados no gráfico aparecem na figura: um gate em outros parâmetros, ou uma combinação, não tem lugar nela. Com várias amostras sobrepostas, os gates são os da primeira amostra, como na tela.

## Cores

As categorias (amostras, grupos, clusters) recebem as cores qualitativas de Paul Tol, que continuam distinguíveis com as deficiências comuns de visão de cores. A ordem das categorias é a ordem escolhida para as amostras e os grupos (veja [[sample-order]]). As escalas para valores são as do matplotlib (viridis, magma, cividis, turbo…), ou uma escala com as suas próprias cores.
