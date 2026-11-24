---
title: Concatenação balanceada
---

O UMAP, o tSNE e os métodos de agrupamento rodam sobre um conjunto de eventos tirados de todas as
amostras: os eventos concatenados. Quantos cada amostra cede decide do que o mapa é mapa. Uma
amostra tratada muitas vezes tem várias vezes as células do seu controle, e um embedding feito
com tudo é um embedding do grupo tratado. No MIL-C o plano — o que cada amostra tem, o que vai
ceder e qual delas limitou — fica na mesa antes de qualquer sorteio, para que o custo do
balanceamento seja lido antes de ser pago.

## Os ajustes

Na seção **Downsampling** da [[samples-tab]]:

| Controle | Significado |
|---|---|
| **Balance** | uma de três maneiras de decidir o número de cada amostra (abaixo) |
| **N** | o número pedido; `as many as the data allows` (0) deixa os dados decidirem |
| **Minimum** | uma amostra com menos eventos da população que isso fica de fora, em vez de puxar todas as outras para o seu tamanho |
| **Seed** | a semente do sorteio; a mesma semente dá os mesmos eventos |

Os eventos são os da **Population** escolhida acima, nas amostras marcadas, agrupadas pelas
colunas de **Group by**.

## As três maneiras de balancear

**Same number from every sample** (o mesmo número de cada amostra). Toda amostra cede o mesmo
número: o da menor amostra, ou **N** se for menor. A Note da menor amostra diz `sets the limit`.
Se **N** for maior do que a menor amostra tem, a nota `asked for N per sample; amostra has only
k` avisa, e a menor continua definindo o número.

**Same total from every group** (o mesmo total de cada grupo). Todo grupo cede o mesmo total: o do
menor grupo, ou **N** se for menor (aqui **N** é por grupo). Dentro de um grupo o total é
repartido tão igualmente quanto as amostras permitem: uma amostra pequena cede tudo o que tem, e
as maiores dividem o que falta, com os eventos que sobram da divisão indo para as maiores. Quando
amostras de um grupo acabam cedendo números diferentes, a nota diz `inside a group the samples do
not give the same number: the larger ones make up for the smaller`.

**At most N from each sample** (no máximo N de cada amostra). Cada amostra cede **N**, ou tudo o
que tem se for menos. Este modo precisa de um número em **N**. Ele não é balanceado quando algumas
amostras ficam abaixo, e avisa: `k sample(s) have fewer than N and give all they have: the
samples are not balanced`.

## Um exemplo

Três amostras Control e três Treated, com estes eventos na população:

| Amostra | Grupo | Tem | Mesmo número | Mesmo total por grupo | No máximo 20.000 |
|---|---|---|---|---|---|
| sample_01 | Control | 12.000 | 12.000 | 12.000 | 12.000 |
| sample_02 | Control | 15.000 | 12.000 | 15.000 | 15.000 |
| sample_03 | Control | 20.000 | 12.000 | 20.000 | 20.000 |
| sample_04 | Treated | 40.000 | 12.000 | 15.666 | 20.000 |
| sample_05 | Treated | 55.000 | 12.000 | 15.667 | 20.000 |
| sample_06 | Treated | 60.000 | 12.000 | 15.667 | 20.000 |
| **total** | | 202.000 | 72.000 | 94.000 | 107.000 |

- Mesmo número: limitado por sample_01; 36.000 por grupo; 36% do que as amostras têm.
- Mesmo total: limitado pelo grupo Control (47.000); cada amostra Treated cede cerca de um terço
  disso.
- No máximo 20.000: os grupos cedem 47.000 e 60.000, e a nota avisa que duas amostras cedem tudo
  o que têm.

Qual usar: **Same number** quando toda amostra deve pesar o mesmo — cada animal é uma observação.
**Same total** quando os grupos devem pesar o mesmo e as amostras são desiguais dentro deles.
**At most N** para limitar arquivos muito grandes sem descartar os pequenos, aceitando que o
resultado não é balanceado. **Minimum** impede que uma amostra quase vazia derrube todas as
outras: com **Minimum** 5.000, uma amostra de 800 eventos fica de fora e a Note dela diz por quê.

## O que é o sorteio

Os eventos de cada amostra são sorteados, sem reposição, dentro da população, e mantidos na
ordem em que foram adquiridos. O sorteio depende da semente e do próprio nome da amostra — não de
quais outras amostras entram nem em que ordem — então desmarcar uma amostra não muda os eventos
tirados das outras.

## O que Concatenate monta

O conjunto concatenado guarda, para cada evento sorteado:

- todos os canais dos arquivos, com a fluorescência compensada (pela matriz que cada amostra
  recebe), nada transformado: a escala é escolhida onde os valores são usados;
- de que amostra ele veio e qual evento daquele arquivo ele era;
- a anotação da amostra (todas as colunas da tabela).

Junto fica registrado como ele foi montado: a população, o modo, o que cada amostra tinha e
cedeu, quais ficaram de fora e por quê, qual amostra ou grupo limitou, e a semente. Onde houve
controle de qualidade ([[quality-control]]) o registro diz quanto foi tirado de cada amostra; onde
a normalização de lotes está ligada ([[normalisation]]) os lotes são alinhados antes de o
conjunto ficar pronto, e o registro diz quanto. Tudo isso aparece no relatório de [[the-plate]] e
é salvo com o conjunto ([[projects-and-saving]]).

Depois de **Concatenate** a janela passa para [[explore-tab]]. Uma nova concatenação substitui o
conjunto anterior e os resultados calculados sobre ele: salve o conjunto antes se quiser
guardá-lo.

## Quando fica desatualizado

Se os valores por baixo das amostras mudarem depois que o conjunto foi montado — outra matriz de
compensação, unmixing, normalização aplicada às próprias amostras — uma linha acima das abas diz
que Explore e Clusters ainda mostram os valores antigos até você concatenar de novo. Mudar o
próprio plano (marcações, **N**, a população) não mexe num conjunto já montado.

## Pela linha de comando

`milc concat` faz o mesmo sem a janela, imprime o plano e grava `NOME.npz` e `NOME.fcs`;
`--dry-run` só imprime o plano:

```
milc concat DADOS --wsp analysis.wsp --population Cells/Live/CD45 \
    --mode equal_per_group --group-by Group --min-events 1000 --seed 0 --out concat
```

Os modos são `equal_per_sample`, `equal_per_group` e `ceiling`. Veja [[command-line]].
