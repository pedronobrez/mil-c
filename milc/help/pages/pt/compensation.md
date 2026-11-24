---
title: Compensação
---

Um fluorocromo é visto não só no seu próprio detector mas, um pouco, nos outros. A matriz de spillover diz quanto: a linha *i* é o fluorocromo *i*, e o valor na coluna *j* é a fração dele que o detector *j* vê. O que o citômetro mede é o sinal verdadeiro vezes a matriz; a compensação multiplica pela inversa da matriz para recuperar o sinal verdadeiro. Nos gates, um parâmetro compensado leva o nome `Comp-` seguido do detector, como o FlowJo escreve.

## De onde vem a matriz de uma amostra

Cada amostra é compensada com uma matriz. Em ordem de precedência:

| Origem | Aparece como | Significado |
|---|---|---|
| Editada aqui para esta amostra | `edited` | Uma matriz mudada no diálogo abaixo e dada só a esta amostra. |
| Editada aqui para todas | `shared` | Uma matriz mudada no diálogo e dada a todas. |
| O workspace do FlowJo | `workspace` | A matriz que o workspace dá a esta amostra, ou a que todas as amostras do workspace compartilham (veja [[workspaces-and-gating-ml]]). |
| O arquivo FCS | `file` | A keyword `$SPILLOVER` (ou `SPILL`) que o citômetro gravou. |
| Nenhuma | `none` | O arquivo não traz matriz: os valores são usados como adquiridos. |

A legenda sob a barra da [[gates-tab]] diz qual delas a amostra mostrada usa. Com o unmixing espectral ligado, nenhuma matriz é usada (veja [[spectral-unmixing]]).

## O diálogo Compensation

**Compensation…** na aba Gates abre o diálogo. Se nenhuma amostra tem matriz, a barra de status diz `These files carry no compensation matrix` e nada abre.

À esquerda fica a lista de amostras, com a origem da matriz de cada uma na coluna **Matrix** (em negrito quando editada aqui). Abaixo, uma linha diz quantas matrizes diferentes as amostras têm, por exemplo `3 different matrices among the 12 samples.`

À direita fica a matriz da amostra selecionada, em porcentagem: cada célula é a porcentagem do fluorocromo da linha vista no detector da coluna. As linhas trazem o detector e o seu marcador. A diagonal é 100 e não pode ser editada. Digite sobre qualquer outro número para mudá-lo; a vírgula é aceita como separador decimal.

Acima da tabela, uma linha diz o que está sendo mostrado: a matriz da própria amostra e de onde ela vem, ou `… matrix with your changes — not applied yet`.

| Botão | O que faz |
|---|---|
| **Apply to** *amostra* / **Apply to the N selected samples** | Dá a matriz da tabela às amostras selecionadas na lista; as outras ficam com a sua. Selecione várias com Shift ou Cmd/Ctrl. |
| **Apply to every sample** | Uma matriz para todas as amostras, no lugar da de cada uma. |
| **Back to the original** | As amostras selecionadas voltam à matriz do seu arquivo, ou do workspace. Fica ativo quando uma delas tem matriz editada. |
| **Discard changes** | Esquece o que foi digitado e mostra de novo a matriz da amostra. |
| **From controls…** | Calcula uma matriz a partir de controles de uma cor e a põe na tabela, para aplicar como uma digitada à mão. Veja [[autocompensation]]. |
| **Close** | Fecha o diálogo. O que foi aplicado continua aplicado. |

Escolher outra amostra na lista enquanto há mudanças ainda não aplicadas mantém as mudanças na tabela, para que a mesma edição possa ir para outras amostras.

Uma matriz é recusada quando uma célula não é número (`Every cell has to be a number`), quando não pode ser invertida (`two fluorochromes would be indistinguishable`), ou quando falta a uma amostra que a recebe um dos seus detectores.

## O que acontece quando uma matriz muda

- Todo valor compensado com ela muda, então todas as contagens são refeitas.
- Os gates foram desenhados nos valores antigos: uma barra acima do gráfico avisa até que se aperte **Reviewed**. Retângulos e polígonos ficam onde estão nos valores novos e podem não estar mais sobre a sua população. Os gates ajustados por **Auto…** são ajustados de novo (veja [[auto-gates]]).
- Os eventos concatenados, se houver, foram tirados dos valores antigos; a janela avisa que estão desatualizados até serem feitos de novo (veja [[concatenation]]).
- As matrizes editadas são salvas com o projeto (veja [[projects-and-saving]]) e vão para um workspace do FlowJo exportado por **File ▸ Export gates…**.

Um arquivo FCS exportado com **Export…** e **Compensated values** desmarcado leva a matriz em uso nas suas keywords; marcado, os valores são gravados já compensados, em colunas chamadas `Comp-…` (veja [[tables-and-prism]]).

## Quando editar

Edite uma matriz quando um gráfico de dois parâmetros de fluorescência mostrar uma população inclinada ao longo de uma diagonal quando deveria estar de pé (subcompensada) ou se curvando para o outro lado (sobrecompensada). **Matrix…** na aba Gates mostra todos os pares de uma vez, que é onde o spillover melhor aparece. Uma mudança digitada de um ou dois pontos percentuais é comum; uma mudança grande em geral quer dizer que os controles deveriam ser adquiridos de novo. Confira também a escala: num eixo biexponencial, um **Linear width** estreito pode fazer uma população negativa bem compensada parecer torta (veja [[scales-and-transforms]]).

As matrizes dos arquivos de dias de aquisição diferentes muitas vezes são diferentes. A linha de resumo do diálogo mostra isso de cara; decidir qual está certa é questão dos controles de cada dia.
