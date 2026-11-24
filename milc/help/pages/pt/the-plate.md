---
title: A prancha
---

Uma prancha é uma página de vários painéis, como uma figura de artigo: um mapa, um heatmap, dois gráficos de gates e um gráfico de frequências, com as letras A, B, C… e gravados num só arquivo. Também pode ser gravada como relatório: um PDF cuja primeira página diz o que foi a análise, seguida das pranchas, uma por página. Um projeto pode ter quantas pranchas forem necessárias, cada uma com um nome.

## Pôr painéis numa prancha

Os painéis vêm do diálogo Figure (veja [[figures]]). Deixe a figura como deve ficar e aperte **Add to plate**, com a prancha escolhida na lista ao lado — uma que já existe, ou **New plate**. A barra de status diz que letra o painel recebeu. Abra a prancha em **View ▸ Plate…**.

Um painel guarda o que mostrava e o acabamento que tinha: rótulos, legenda, cores e escala de cores. A página decide a fonte, o tamanho e a resolução de todos os painéis, que é o que faz de seis gráficos uma figura só.

## A janela Plate

A prévia à esquerda é desenhada pelo mesmo código que grava o arquivo. À direita:

| Parte | O que faz |
|---|---|
| **Plate** | A lista das pranchas do projeto, e **New**, **Rename…**, **Duplicate** e **Delete**. Uma prancha nova tem a fonte e o tamanho da atual e nenhum painel. **Duplicate** faz uma cópia para mudar separadamente da original. **Delete** pergunta antes. |
| **Panels, in order** | Os painéis da prancha, com as letras. Um painel cujo gráfico não tinha nada a mostrar no último refresh aparece marcado como desatualizado (stale). |
| **Move up**, **Move down**, **Remove** | Reordenam os painéis, ou tiram um. |
| **Refresh from the current analysis** | Desenha cada painel de novo a partir do gráfico de onde veio, como esse gráfico está agora. Um painel cujo gráfico não tem nada a mostrar fica com o que tinha e é marcado como desatualizado. |

**The panel selected** (o painel selecionado):

| Controle | O que faz |
|---|---|
| **Title** | Um título sobre o painel; nenhum por padrão. |
| **Takes** | Quantas colunas e linhas da grade ele ocupa, de 1 a 6 cada. |

**Page** (a página):

| Controle | Padrão | O que faz |
|---|---|---|
| **Panels per row** | 2 | Colunas da grade, de 1 a 6. |
| **Fit** | | **Single column, 85 mm**, **1.5 columns, 120 mm**, **Double column, 180 mm**, ou **Custom**. |
| **Width × height** | 180 mm × Automatic | A altura **Automatic** é três quartos da largura de uma célula para cada linha. |
| **Letters: A, B, C…** | ligado | Uma letra em negrito no canto superior esquerdo de cada painel, 4 pontos maior que o texto. Depois de Z vêm AA, AB… |
| **Font**, **Size**, **Resolution** | | A fonte de todos os painéis, e a resolução de um PNG ou TIFF. |

Os painéis são postos em ordem, cada um no primeiro lugar livre que o comporta, de modo que as letras se leiam da esquerda para a direita e depois para baixo, seja qual for o tamanho dos painéis.

## Gravar

- **Export…** grava a prancha atual como **SVG**, **PDF**, **PNG** ou **TIFF**, pela extensão do arquivo.
- **Report…** grava um PDF. Se várias pranchas têm painéis, primeiro pergunta quais entram, uma por página. Quando há eventos concatenados, a primeira página (ou as primeiras) diz o que foi a análise, a partir do que os dados concatenados registram: quantos eventos de quantas amostras, a população, o controle de qualidade, a amostragem (modo, semente, quanto se tirou de cada amostra e o que ficou de fora), a normalização por lote e as execuções (UMAP, tSNE, FlowSOM e outras) com os seus ajustes. Nada ali é recalculado nem interpretado. O PDF registra a versão do MIL-C que o gravou.

## Salva com o projeto

As pranchas são gravadas com o projeto como estão, no `MIL-C.json` da pasta e, para os seus arrays, no `MIL-C plates.npz` ao lado (veja [[projects-and-saving]]). Um painel lido de volta de um projeto é o que era quando foi salvo: não segue mais um gráfico, e **Refresh from the current analysis** o deixa como está. É de propósito: a página de um artigo não deve mudar porque a análise por trás dela rodou mais uma vez. Para atualizar um painel desses, faça a figura de novo e acrescente-a.

Um painel que não pode ser lido de volta (um tipo de gráfico que esta versão não conhece) aparece como ausente na página, e o que estava guardado é mantido, então salvar de novo não perde nada.
