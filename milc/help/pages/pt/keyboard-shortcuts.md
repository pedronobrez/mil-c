---
title: Atalhos de teclado
---

O MIL-C é usado sobretudo com o mouse: gates são desenhados, amostras clicadas, valores
arrastados. O teclado tem poucos atalhos, listados aqui junto com os gestos de mouse que não são
óbvios pelos botões. Os atalhos seguem as convenções de cada sistema, então `Ctrl` no Windows e no
Linux é `⌘` no macOS.

## A janela principal

| Ação | Windows, Linux | macOS | Menu |
|---|---|---|---|
| Abrir uma pasta de `.fcs` | `Ctrl+O` | `⌘O` | `File ▸ Open folder…` |
| Salvar os eventos concatenados | `Ctrl+S` | `⌘S` | `File ▸ Save concatenated file…` |
| O manual, na página de onde você está | `F1` | `⌘?` | `Help ▸ Manual` |

O `F1` funciona também dentro de um diálogo: no diálogo de qualidade ele abre
[[quality-control]], no de figura [[figures]], no de compensação [[compensation]]. Onde um
diálogo não tem página própria, o manual abre na página da aba por trás dele.

Os outros comandos dos menus — `File ▸ Open workspace…`, `Open concatenated file…`, `Export
tables…`, `Export gates…`, `View ▸ Appearance`, `View ▸ Plate…`, `Help ▸ Export manual as
PDF…` — não têm atalho.

## A janela do manual

| Ação | Windows, Linux | macOS |
|---|---|---|
| Buscar no manual | `Ctrl+F` | `⌘F` |
| Abrir o primeiro resultado | `Return` na caixa de busca | `Return` |
| Voltar à página anterior | `Alt+←` | `⌘[` |
| Avançar de novo | `Alt+→` | `⌘]` |

Veja [[welcome]] para como a busca casa as palavras.

## A tabela da aba Samples

| Ação | Como |
|---|---|
| Selecionar várias linhas | `Shift`-clique para uma sequência, `Ctrl`-clique (`⌘`-clique) para mais uma |
| Editar uma célula das suas colunas | duplo clique, ou selecione e comece a digitar; `F2` no Windows e no Linux |
| Pôr um valor em muitas amostras | selecione as linhas delas e digite numa |
| Cancelar uma edição | `Esc` |
| Renomear uma coluna sua | duplo clique no cabeçalho |
| Renomear ou remover uma coluna | botão direito no cabeçalho |

Veja [[samples-tab]].

## Desenhando gates

No gráfico da [[gates-tab]], com uma ferramenta de gate escolhida:

| Ferramenta | Gesto |
|---|---|
| **Rectangle** | arraste sobre os eventos; num histograma, uma faixa |
| **Polygon** | clique em cada vértice; duplo clique para fechar |
| **Lasso** | segure o botão do mouse e desenhe em volta dos eventos; solte para fechar |
| **Ellipse** | arraste a caixa em que ela cabe, depois gire e estique pelas alças |
| **Quadrant** | clique onde as duas linhas devem se cruzar |
| qualquer uma | `Esc` interrompe o desenho |

A lista de amostras da aba Gates aceita `Shift`-clique e `Ctrl`-clique (`⌘`-clique) para
selecionar várias amostras e sobrepô-las. Mais em [[drawing-gates]].

## Mapas e diálogos

- No mapa da [[clusters-tab]], `Shift`+arrastar contorna à mão livre qualquer conjunto de eventos,
  e a aba dá para ele os mesmos números que dá para um cluster.
- No diálogo de figura, duplo clique na cor de uma categoria para mudá-la; digite por cima do
  rótulo para renomeá-la.
- No diálogo Order of the groups, arraste um valor para movê-lo ([[sample-order]]).
- Os diálogos fecham com `Esc` e confirmam com `Return` quando têm um botão **OK**, como em todo o
  sistema.
