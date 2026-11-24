---
title: A aba Gates
---

A aba Gates mostra uma amostra por vez: dois parâmetros quaisquer dela, ou o histograma de um, para uma população da árvore de gates. É onde os gates são desenhados, movidos, renomeados e apagados, onde se escolhem a compensação e a escala de cada eixo, e onde se ajustam os modelos de um parâmetro (proliferação, ciclo celular, cinética).

Um gate desenhado aqui pertence à população que está sendo vista, como no FlowJo. Com **Apply to every sample** marcado, ele pertence também a todas as outras amostras, do jeito que um gate posto num grupo do FlowJo. O que se desenha fica guardado com o projeto, ao lado dos arquivos (veja [[projects-and-saving]]), e passa a valer no lugar dos gates do workspace.

## O lado esquerdo

| Parte | O que é |
|---|---|
| **Samples** | Os arquivos utilizáveis da pasta, cada um com seu número de eventos. Clique em um para vê-lo. Selecione vários (Shift-clique, ou Cmd/Ctrl-clique) para sobrepô-los. |
| **Populations** | A árvore de gates da amostra mostrada, com três colunas: **Population**, **Events** e **% parent**. Clique numa população para vê-la. |
| **New gate on the population shown** | Os nove botões que desenham e editam gates, descritos em [[drawing-gates]]. |
| **Apply to every sample** | Marcado (o padrão), um gate desenhado, movido, renomeado ou apagado nesta amostra é desenhado, movido, renomeado ou apagado em todas. Desmarcado, só nesta. |
| **Back to the workspace's gates** | Esquece todo gate desenhado ou movido aqui e volta aos gates do workspace do FlowJo, ou a nenhum gate. Fica ativo depois que algo foi mudado. |

Uma população em *itálico* foi lida de um workspace do FlowJo por uma regra que ainda não foi conferida com o próprio FlowJo; a dica dela mostra a contagem do FlowJo ao lado da contagem daqui. Veja [[workspaces-and-gating-ml]].

Ao escolher uma população, os eixos passam para os dois parâmetros do primeiro gate dela, de modo que a população é vista primeiro do jeito que seus gates foram desenhados. Depois de desenhar um gate, os eixos ficam onde estão.

## A barra acima do gráfico

| Controle | O que faz |
|---|---|
| **X**, **Y** | Os parâmetros de cada eixo, com o nome do marcador quando o arquivo traz um. A lista de **Y** termina em **Histogram**, que mostra só o parâmetro de X. |
| **Scale…** (ao lado de X e de Y) | A escala daquele eixo: biexponencial, Logicle, arcsinh, logarítmica ou linear. Veja [[scales-and-transforms]]. |
| **Plot** | **Pseudocolour** ou **Contour**, para dois parâmetros e uma amostra. |
| **Colour by** | Com várias amostras sobrepostas: uma cor para cada **Sample**, ou para cada valor de uma das suas colunas (veja [[samples-tab]]). |
| Escolha da sobreposição | Com várias amostras num histograma: **Overlaid**, **Staggered**, **Violins** ou **Boxes**, e **One for each group** para juntar as amostras de cada valor de **Colour by** num só violino ou caixa. |
| **Quality…** | Procura os trechos da aquisição que deram errado. Veja [[quality-control]]. |
| **Compensation…** | A matriz de spillover de cada amostra. Veja [[compensation]]. |
| **Unmixing…** | Unmixing espectral no lugar da compensação. Veja [[spectral-unmixing]]. |
| **Export…** | Os eventos da população mostrada, em arquivos FCS 3.1 ou CSV. Veja [[tables-and-prism]]. |
| **Figure…** | O gráfico como figura para artigo. Veja [[figures]]. |

Abaixo da barra, uma legenda diz o que está desenhado, por exemplo `sample_01 · CD4+: 12,034 events, 41.2% of Lymphocytes — 40,000 drawn · compensated with its file's matrix`. Ela também avisa quando há amostras sobrepostas (e de quem são os gates mostrados: os da primeira amostra), quando os valores estão normalizados por lote e qual matriz compensa a amostra.

## O que o gráfico mostra

- **Pseudocolour**: cada evento é um ponto colorido pela densidade da vizinhança, numa escala logarítmica (uma grade de 160 × 160 sobre o gráfico). Os pontos mais densos ficam por cima. Acima de 40.000 eventos, desenham-se 40.000 sorteados (sempre os mesmos); todas as contagens e porcentagens continuam usando todos os eventos.
- **Contour**: linhas que envolvem 90, 75, 50, 25 e 10% dos eventos; os eventos fora da linha mais externa aparecem como pontos.
- **Histogram**: 256 classes ao longo do eixo, com a contagem de eventos à esquerda. Com várias amostras, cada uma é escalada pela sua classe mais alta (`% of the tallest bin`); **Staggered** põe uma acima da outra.
- **Violins** e **Boxes**: a distribuição do parâmetro de X em cada amostra, uma linha para cada, com o eixo na horizontal como no histograma, para que um gate de faixa fique no mesmo lugar.

Cada gate da população mostrada aparece no gráfico com o nome e a fração da população, quando os dois parâmetros do gráfico são os do gate. Gates em outros parâmetros não aparecem no gráfico, mas continuam na árvore.

## Ao lado da legenda

| Botão | O que abre |
|---|---|
| **Matrix…** | Todos os pares de parâmetros da população mostrada, numa grade só: é onde o spillover aparece. Marque os parâmetros na lista (no máximo 12; desenham-se os primeiros marcados). Começa com os parâmetros de fluorescência marcados. Cada célula é colorida pela densidade de eventos, numa escala logarítmica própria. |
| **Backgate…** | A população mostrada em cada gate acima dela: em cinza, todos os eventos aos quais o gate foi aplicado; em cor, os eventos da população, onde quer que caiam. Não se aplica a **All events**, que não tem gate acima. |
| **Model…** | Um menu: **Proliferation…**, **Cell cycle…** e **Kinetics…**, ajustados à população mostrada em cada amostra selecionada. Veja [[proliferation]], [[cell-cycle]] e [[kinetics]]. |

**Matrix…** e **Backgate…** têm cada um o seu botão **Figure…**.

## Avisos acima do gráfico

Duas barras aparecem sobre o gráfico quando os gates podem não estar mais onde deveriam.

- **These gates were drawn before …** aparece depois que algo mudou os valores que todo gate corta: uma nova matriz de compensação, o unmixing, ou a normalização por lote aplicada aos gates. Olhe os gates, mova o que precisar, e aperte **Reviewed**.
- **… was fitted automatically and then moved by hand …** aparece quando um gate proposto por **Auto…** foi movido à mão e a população em que ele está mudou desde então. **Refit** ajusta de novo na população como está agora e abre mão do ajuste manual; **Keep** deixa onde foi posto e para de avisar. Veja [[auto-gates]].

## Teclas

**Esc** cancela o gate que está sendo desenhado. Os outros atalhos da janela estão em [[keyboard-shortcuts]].
