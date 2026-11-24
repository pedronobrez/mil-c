---
title: A aba Clusters
---

A aba Clusters mostra os clusters de quatro lados ao mesmo tempo: o **Map** diz onde um cluster está, o heatmap **Expression** o que ele expressa, as barras de **Frequency** em quais amostras e grupos ele aparece, e a **FlowSOM tree** como os nós do mapa se relacionam. Um clique num cluster em qualquer deles o marca nos quatro. Ela fica vazia até o FlowSOM ou o PhenoGraph serem rodados na [[explore-tab|aba Explore]]; como eles funcionam está em [[clustering-methods]].

## A primeira linha: a aparência

| Controle | O que faz |
|---|---|
| **Clusters** | Qual clusterização é mostrada, FlowSOM ou PhenoGraph. Só aparece quando as duas foram rodadas |
| **Map** | Sobre o que o mapa é desenhado: UMAP ou tSNE, onde foram rodados, ou os dois primeiros marcadores em que a clusterização rodou (por exemplo `CD3 × CD4`) |
| **Heatmap** | **Median** ou **Mean** de cada marcador em cada cluster |
| **Colours** | A rampa de cores do heatmap, da árvore colorida por um marcador e da expressão nos mapas: uma rampa com nome ou duas ou três cores suas |
| **Colour per marker** | Marcado (padrão): cada coluna do heatmap percorre a rampa inteira por si, para que um marcador fraco mostre as suas próprias diferenças. Desmarcado: uma rampa para a tabela toda |
| **Numbers** | O valor de cada célula do heatmap escrito nela, aqui e na figura |
| **Other markers too (n)** | Mostra também os marcadores corados em que a clusterização não rodou. Diz "(none)" e fica desligado quando ela rodou em todos |
| **Tree** | Colore a árvore do FlowSOM por **Metacluster** ou por qualquer marcador |
| **Clear selection** | Solta o cluster ou o laço selecionado |

## A segunda linha: o que se compara

| Controle | O que faz |
|---|---|
| **Compare** | A coluna das amostras cujos grupos as barras comparam por cor. Começa no agrupamento escolhido em Samples. Com várias colunas escolhidas lá, a combinação delas também é oferecida. "All samples" quando as amostras não têm coluna que varie |
| **within each** | Uma segunda subdivisão: um painel de barras para cada um dos seus valores, de modo que os grupos são comparados dentro de cada um — os grupos dentro de cada dia, por exemplo |
| **Order…** | Qual grupo vem primeiro, em todo gráfico, legenda, painel e tabela: ver [[sample-order]] |
| **Test** | **No test**, ou um teste da fração dos grupos em cada cluster: **Welch t test / ANOVA**, **Mann-Whitney / Kruskal-Wallis**, **Welch t test / Welch ANOVA** |
| (correção) | Ao lado do teste: **Benjamini-Hochberg**, **Holm**, **Bonferroni** ou **None** |
| **Paired by** | A coluna que diz de qual animal ou doador é uma amostra, quando os mesmos foram medidos em todos os grupos |
| **Results…** | Todas as comparações numa tabela, par a par quando há mais de dois grupos |
| **Prism…** | A fração de cada cluster em cada amostra como arquivo do GraphPad Prism: ver [[tables-and-prism]] |
| **Differential abundance…** | Dois grupos comparados pelas contagens de células dos clusters: ver [[abundance]] |
| **Differential expression…** | Dois grupos comparados pelo que as células de cada cluster expressam: ver [[differential-expression]] |

O que os testes fazem e como lê-los está em [[cluster-statistics]].

## Os quatro painéis

Cada painel tem o seu título e um botão **Figure…** que o transforma em figura (ver [[figures]]). As divisórias entre eles podem ser arrastadas.

- **Map** — cada evento concatenado, colorido pelo seu cluster. Com um cluster selecionado, os outros ficam cinza e continuam onde estão, para que o cluster seja lido contra o mapa inteiro.
- **Expression** — clusters nas linhas, marcadores nas colunas; cada célula é a mediana (ou média) do marcador sobre todos os eventos concatenados do cluster, na escala em que a clusterização rodou (0 a 1023 na biexponencial). A cor vai por marcador, a menos que **Colour per marker** esteja desmarcado.
- **Frequency** — para cada cluster, a fração dos eventos de cada amostra que cai nele: uma barra por grupo na média, uma haste de um desvio-padrão, e cada amostra como um ponto — cinco animais são cinco pontos, não só uma barra. A legenda dá o n de cada grupo. Com **within each**, um painel para cada um dos seus valores, na mesma escala. Com um teste escolhido, o p ajustado é escrito acima de cada cluster onde fica abaixo de 0,05, e "ns" onde não fica.
- **FlowSOM tree** — um nó para cada unidade do mapa que tem eventos, do tamanho da sua população, ligados pela árvore geradora mínima dos nós; colorida por metacluster ou por um marcador. Só aparece para o FlowSOM.

## Selecionar

- **Clique** num cluster no heatmap, nas barras ou na árvore para selecioná-lo; clique de novo para soltá-lo.
- **Shift+arrastar** no mapa contorna qualquer conjunto de eventos (um laço).

A linha abaixo dos painéis descreve a seleção: o número de eventos e a fração de todos, a nota do cluster, e a porcentagem média da amostra em cada grupo (± desvio-padrão, n). Para um laço, ela também diz a quais três clusters pertence a maioria dos seus eventos.

## Trabalhar os clusters

A fileira de botões embaixo age sobre o cluster selecionado. O resultado do algoritmo nunca é sobrescrito: nomes, notas, cores e junções ficam por cima dele, salvos com a execução, de modo que uma junção sempre pode ser desfeita. Uma nova execução do algoritmo começa do zero.

| Botão | O que faz |
|---|---|
| **Rename…** | Dá um nome ao cluster selecionado. O nome aparece em todo lugar onde estava o número. Vazio: volta ao número. Dois clusters não podem ter o mesmo nome |
| **Note…** | Uma nota guardada com o cluster, mostrada na linha abaixo dos painéis |
| **Colour…** | Uma cor à sua escolha para o cluster no mapa, na árvore e nas figuras deles, ou de volta à cor da própria paleta |
| **Merge with…** | Faz um cluster só do selecionado e de outros marcados numa lista. Ele fica com o menor número e mantém o nome e a nota desse |
| **Undo merges** | Volta aos clusters que o algoritmo deu; nomes e notas dos clusters juntados vão embora também |
| **Suggest names…** | Propõe um nome para cada cluster a partir dos marcadores em que ele é positivo — `CD4+ Foxp3+` — para aceitar, editar ou recusar |
| **Make population** | O cluster do FlowSOM selecionado como população dos arquivos, em Gates e Statistics |
| **Save selection…** | Guarda os eventos do laço como um cluster próprio; os outros eventos ficam nos seus clusters |
| **Remove** | Remove o cluster desenhado à mão selecionado: os seus eventos voltam aos clusters em que estavam |

### Como os nomes são sugeridos

Para cada marcador em que a clusterização rodou, os valores dos clusters são ordenados e cortados em dois no maior intervalo entre vizinhos: os de cima são positivos, os de baixo negativos. Um marcador cujo maior intervalo é mais estreito que um décimo da escala biexponencial de 0 a 1023 (cerca de meia década), ou 1 numa escala arcsinh, não separa nada e fica de fora. Um cluster positivo para nenhum recebe o nome dos marcadores em que é negativo (`CD4- CD8-`); dois que os marcadores não distinguem recebem ` (2)`, ` (3)` depois do nome. O diálogo mostra o nome atual de cada cluster ao lado da sugestão; escreva por cima de qualquer um, ou apague um para deixar aquele cluster como está, e depois **Apply**.

### Um cluster como população

**Make population** cria, em cada arquivo, uma população com o nome do cluster (`FlowSOM 3` para um numerado) sob a população de onde a concatenação foi tirada. Cada evento de cada arquivo — não só os sorteados para a concatenação — vai para o seu nó mais próximo do mapa treinado do FlowSOM, e pertence à população se esse nó é um dos do cluster. A população aparece então em [[gates-tab|Gates]] e [[statistics-tab|Statistics]] como qualquer gate. A operação é recusada para um cluster desenhado à mão (o mapa em que foi contornado só tem coordenadas para os eventos sorteados), para clusters do PhenoGraph, quando o mapa treinado não está na memória (um arquivo concatenado aberto do disco não o traz: rode o FlowSOM de novo) e quando o mapa foi treinado com eventos normalizados por lote mas os arquivos não são normalizados.
