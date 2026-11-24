---
title: UMAP, tSNE e PCA
---

Uma redução de dimensão põe cada evento num plano de modo que eventos parecidos em muitos marcadores fiquem perto uns dos outros. A [[explore-tab|aba Explore]] oferece três: UMAP e tSNE, que desenham o mapa que a maioria espera, e PCA, cujos dois eixos podem ser lidos. As três rodam sobre os eventos concatenados (ver [[concatenation]]), nos marcadores marcados e na escala escolhida em Explore.

Os padrões são os dos plugins do FlowJo, de modo que um primeiro mapa é comparável ao que ele substitui.

## O que os algoritmos recebem

| Escala | Valores |
|---|---|
| Biexponential (padrão) | Cada marcador na escala que o seu eixo tem em Gates — a biexponencial do FlowJo, a menos que outra tenha sido escolhida lá — como números de 0 a 1023. É o que o FlowJo entrega aos seus plugins. Um marcador mostrado em escala linear em Gates vai para 0 a 1023 linearmente; um mostrado em arcsinh em Gates vai como o seu arcsinh |
| Arcsinh, cofactor 150 | asinh(x / 150) de todos os marcadores |
| Logicle, W 0.5, M 4.5 | A Logicle padrão (topo 262.144) de todos os marcadores, em 0 a 1023 |

Nada mais é reescalado ou centrado, exceto que o PCA centra os valores. A fluorescência é compensada (ver [[compensation]]). A escala importa: numa escala que achata os negativos, eventos perto de zero se fundem; numa linear, marcadores brilhantes dominam. Ver [[scales-and-transforms]].

## UMAP

O UMAP (McInnes et al., 2018) monta um grafo dos vizinhos mais próximos de cada evento no espaço dos marcadores e o estende no plano de modo a preservar o grafo. O MIL-C chama o pacote `umap-learn`, com distância euclidiana, duas dimensões e o seu layout original ("compatibility").

| Configuração | Padrão | Efeito |
|---|---|---|
| **Neighbours** | 15 | Com quantos vizinhos cada evento é ligado. Menos: estrutura local mais fina, ilhas mais espalhadas. Mais: um quadro mais amplo, em que as populações grandes mantêm suas posições relativas. Nunca mais que o número de eventos menos um |
| **Min. distance** | 0.5 | Quão juntos os eventos podem ficar no mapa. Menor: ilhas apertadas e compactas. Maior: manchas espalhadas por igual |
| **Seed** | 0 | Com semente o UMAP roda numa só thread, então é reprodutível e mais lento |

## tSNE

O tSNE (van der Maaten e Hinton, 2008) reproduz, no plano, a probabilidade de dois eventos serem vizinhos no espaço dos marcadores. O MIL-C chama o pacote `openTSNE`, distância euclidiana, 1.000 iterações (250 de exagero inicial e depois 750), em todos os núcleos.

| Configuração | Padrão | Efeito |
|---|---|---|
| **Perplexity** | 30 | Aproximadamente o número de vizinhos que cada evento considera. Menor: grupos pequenos e apertados; maior: estrutura mais suave e maior. Limitada a um terço do número de eventos |
| **Seed** | 0 | O arranjo inicial aleatório |

As distâncias entre ilhas num mapa de tSNE dizem pouco; os tamanhos das ilhas também.

## PCA

Os componentes principais são as direções em que os eventos mais variam, achadas por decomposição em valores singulares dos valores centrados. Ao contrário de UMAP e tSNE, saem iguais toda vez, e cada eixo pode ser lido: é tanto deste marcador e tanto daquele.

- Os eixos levam no rótulo a sua fração da variância total, `PCA1 (41%)`.
- Abaixo dos botões, cada componente é descrito pelos seus três marcadores de maior peso (loadings): `PCA1, 41%: +0.72 CD4, −0.55 CD8, +0.30 CD25`. Um marcador com peso positivo grande é alto à direita daquele eixo.
- O sinal de um componente é arbitrário; o MIL-C torna positivo o seu marcador mais pesado, para que uma nova execução aponte para o mesmo lado.
- Não há configuração nem semente. São calculados dois componentes.

## Quão bom é um mapa

Nem o UMAP nem o tSNE dão a mesma imagem duas vezes sem a mesma semente, e nunca a mesma imagem que outra implementação. O que se pode comparar é quão bem um mapa mantém os vizinhos que um evento tinha no espaço dos marcadores. Depois de cada UMAP e tSNE, a aba Explore relata:

> UMAP keeps 47% of each event's 15 nearest neighbours

Isto é: para cada evento, dos seus 15 vizinhos mais próximos nos marcadores, a fração que continua entre os seus 15 mais próximos no mapa, em média sobre os eventos (sobre 20.000 sorteados quando há mais). 1 seria um mapa perfeito; o acaso dá cerca de 15 dividido pelo número de eventos. Os mapas do próprio FlowJo dos mesmos eventos mantêm mais ou menos a mesma fração (ver [[measured-facts]]).

## Como ler um mapa

- As ilhas são grupos de eventos parecidos nos marcadores usados; a posição de uma em relação à outra só tem significado aproximado.
- Um mapa montado com tudo é dominado pelos grupos com mais eventos. A amostragem balanceada da [[concatenation|concatenação]] existe para que toda amostra contribua igualmente.
- Colora por um marcador para ver o que separa as ilhas, e por amostra ou grupo para ver se uma ilha é de uma amostra só.
- A clusterização não é feita sobre o mapa: FlowSOM e PhenoGraph rodam sobre os próprios marcadores. Ver [[clustering-methods]].
