---
title: FlowSOM e PhenoGraph
---

O MIL-C tem duas maneiras de separar os eventos concatenados em clusters, as duas rodadas na [[explore-tab|aba Explore]] e mostradas na [[clusters-tab|aba Clusters]]. Ao **FlowSOM** se diz quantos clusters fazer; o **PhenoGraph** acha o número sozinho. Os dois rodam sobre os marcadores marcados em Explore, na escala escolhida lá (ver [[embeddings]] para o que os algoritmos recebem), e não sobre um mapa de UMAP ou tSNE. Nada mais é reescalado: como nos plugins do FlowJo, os marcadores já chegam numa escala única.

Todo evento termina em exatamente um cluster, numerado a partir de 1, o maior primeiro, para que uma nova execução com outra semente não embaralhe os nomes à toa.

## FlowSOM

O FlowSOM (Van Gassen et al., *Cytometry A* 2015) foi escrito seguindo os passos do próprio pacote R.

1. **O mapa auto-organizável.** Uma grade de nós (10 × 10 por padrão) é posta no espaço dos marcadores, começando em eventos sorteados. Os eventos são apresentados um a um, dez passadas por todos; cada um puxa para si o nó mais próximo e os nós vizinhos dele na grade. A vizinhança é uma bolha que encolhe do quantil 0,67 das distâncias da grade até zero, e a taxa de aprendizado cai de 0,05 para 0,01.
2. **Cada evento no seu nó.** Todo evento vai para o nó mais próximo.
3. **Metaclusters.** Os nós são juntados no número de clusters pedido por clusterização hierárquica de consenso, como faz o ConsensusClusterPlus: cem reamostragens de 90% dos nós, cada uma cortada por ligação média; a frequência com que dois nós caem juntos vira a sua semelhança, e essa é cortada mais uma vez, por ligação média, nos clusters finais.

O cluster de cada evento é o metacluster do seu nó.

| Configuração | Padrão | Efeito |
|---|---|---|
| **Clusters** | 5 | O número de metaclusters, de 2 a 60. Peça mais do que espera: juntar dois na aba Clusters é fácil, separar um não é |
| **Grid** | 10 | O mapa tem este número de nós de cada lado (100 nós com 10). Grades maiores veem estrutura mais fina e levam mais tempo |
| **Seed** | 0 | Os nós iniciais, a ordem dos eventos e as reamostragens |

**A árvore.** A aba Clusters desenha a árvore geradora mínima dos nós que têm eventos, cada nó do tamanho da sua população. O arranjo põe dois nós tão afastados quanto o número de passos entre eles ao longo da árvore — o que o Kamada-Kawai faz no pacote R — e sai igual toda vez.

Os números aleatórios do R não são reproduzidos, então os clusters coincidem com os do FlowJo em substância e não evento a evento: num experimento real, o índice de Rand ajustado entre os clusters do MIL-C e os do FlowJo foi 0,86 e 0,81, contra 0,81 e 0,79 entre duas execuções do MIL-C com sementes diferentes (ver [[measured-facts]]).

Só um cluster do FlowSOM pode virar população dos arquivos completos (**Make population** na aba Clusters), porque o mapa treinado sabe situar qualquer evento.

## PhenoGraph

O PhenoGraph (Levine et al., *Cell* 2015) acha os clusters como comunidades de um grafo. Foi escrito a partir dos artigos, não do pacote PhenoGraph: os passos são os dele, os rótulos não prometem ser os dele.

1. **Vizinhos.** Cada evento é ligado aos seus *k* vizinhos mais próximos por distância euclidiana: de forma exata (uma árvore k-d) até 50.000 eventos ou 10 marcadores, por descida de vizinhos mais próximos além disso, numa só thread para que a mesma semente dê os mesmos vizinhos.
2. **Pesos de Jaccard.** Cada ligação pesa a fração das vizinhanças dos dois eventos que eles têm em comum. Um evento conta na sua própria vizinhança, e um par é ligado quando qualquer um dos dois lista o outro, com o peso inteiro.
3. **Comunidades.** O grafo é cortado nas comunidades que maximizam a modularidade, pelo algoritmo de **Leiden** (Traag et al., 2019) ou pelo método de **Louvain** (Blondel et al., 2008), uma vez, a partir da semente.
4. **Clusters pequenos.** Um cluster com menos eventos que o mínimo é absorvido pelo cluster a que está mais ligado, os menores primeiro. O pacote original deixa esses eventos sem cluster; aqui todo evento precisa estar num cluster, por causa das tabelas de frequência.

| Configuração | Padrão | Efeito |
|---|---|---|
| **Neighbours** (*k*) | 30 | Cada evento é ligado a este número. Menos acha clusters menores; mais os junta |
| **Resolution** | 1.0 | Acima de 1, mais clusters e menores; abaixo de 1, menos e maiores |
| **Method** | Leiden | O Leiden acha um corte um pouco melhor e nunca deixa um cluster em pedaços desconectados; o Louvain é o que o PhenoGraph original usava |
| **Min. cluster size** | 10 | Um cluster com menos eventos é juntado ao vizinho; 1 mantém todos |
| **Seed** | 0 | A busca de vizinhos acima de 50.000 eventos e a ordem em que as comunidades são procuradas |

Depois de rodar, a aba Explore relata o número de clusters e a modularidade, com o modo como foi feito, por exemplo `PhenoGraph: 14 clusters, modularity 0.85 (Leiden, k 30, resolution 1, 2 of fewer than 10 events merged)`.

### Por que pode dar mais clusters do que você espera

A modularidade tem uma escala própria e nenhum modo de dizer "uma população". Com resolução 1, uma população grande sai em vários pedaços, tenha ela partes ou não: uma gaussiana de 3.500 eventos em 8 marcadores sai em 13 clusters com k = 30, aqui e no pacote original igualmente. Quanto melhor o otimizador, mais certamente a população grande é cortada. A dica abaixo das configurações diz o que fazer:

- baixar **Resolution** (0.5 é um bom começo), com o risco de perder uma população rara num vizinho;
- ou aumentar **Neighbours** (100 a 150), a um custo maior de tempo;
- ou manter os pedaços e juntá-los na aba Clusters.

Uma população muito menor que *k* é absorvida pelo vizinho com qualquer resolução.

## Qual usar

| | FlowSOM | PhenoGraph |
|---|---|---|
| Número de clusters | Você escolhe | Achado pelos dados e pela resolução |
| Tempo | Cresce com eventos × nós × passadas | Dominado pela busca de vizinhos, que cresce com eventos e marcadores |
| Reprodutibilidade | Mesma semente, mesmo resultado | Mesma semente, mesmo resultado |
| Como o plugin do FlowJo | Sim, em substância | Sem equivalente no FlowJo |
| Pode virar população dos arquivos | Sim | Não |

Os dois podem ser rodados nos mesmos eventos; a aba Clusters oferece então a escolha entre eles. Nenhum é a verdade: duas execuções de qualquer um, com sementes diferentes, concordam entre si mais ou menos tanto quanto cada um concorda com o pacote original (ver [[measured-facts]]). Olhe as [[cluster-statistics|estatísticas dos clusters]] com sementes diferentes antes de confiar num cluster pequeno.
