---
title: A linha de comando
---

Tudo o que a janela faz com uma pasta pode ser feito sem ela: num script, num servidor, ou numa
execução que precisa ser repetida exatamente. A linha de comando usa a mesma sessão que a janela —
os mesmos gates, matrizes, arquivo de projeto e controle de qualidade — então os números são os
mesmos da janela.

## Como chamar

| Instalado como | Comando |
|---|---|
| código-fonte, ou `pip install` | `milc …` ou `python -m milc …` |
| `.deb` ou `install.sh` no Linux | `milc …` |
| Windows | `MIL-C.exe …`, na pasta do programa |
| macOS | `/Applications/MIL-C.app/Contents/MacOS/MIL-C …` |

Sem comando nenhum, a janela abre. `milc --version` imprime a versão; `milc --help` lista os
comandos, e `milc COMANDO --help` as opções de um deles. Abaixo, `DADOS` é uma pasta de arquivos
`.fcs`.

Um comando que não consegue fazer o que foi pedido imprime o motivo na saída de erro e termina com
status 2.

## inspect: o que uma pasta contém

```
milc inspect DADOS [-v]
```

Uma linha por arquivo: nome, eventos, parâmetros, citômetro e data, só pelos cabeçalhos.
`-v` / `--verbose` lista também cada parâmetro com o seu marcador.

## counts: cada gate contra o FlowJo

```
milc counts DADOS --wsp analysis.wsp [--listed-only] [--pattern REGEX] [--qc …]
```

A árvore de gates de cada amostra com a contagem do MIL-C de cada população ao lado da do FlowJo,
a diferença, e no fim a maior diferença entre todos os gates. Os arquivos que o workspace não cita
são nomeados no fim.

Opções comuns a `counts` e `concat`:

| Opção | Significado |
|---|---|
| `--wsp ARQUIVO` | o workspace do FlowJo (obrigatório) |
| `--listed-only` | só os arquivos que o workspace cita; por padrão os outros recebem os gates e a matriz que estes compartilham |
| `--pattern REGEX` | colunas lidas do nome de cada arquivo; cada grupo nomeado é uma coluna: `'(?P<Day>\d+)(?P<Group>[A-Za-z]+)'`. Sem ela, uma coluna, Group, adivinhada dos nomes ou como digitada na janela |

## concat: subamostragem balanceada e concatenação

```
milc concat DADOS --wsp analysis.wsp [--population CAMINHO] [--mode MODO] [--n N]
     [--min-events N] [--group-by COLUNAS] [--seed N] [--dry-run] [--out NOME]
     [--normalize-by COLUNA --reference COLUNA=VALOR [--normalize-clusters N]
      [--normalize-to LOTE]] [--qc …]
```

| Opção | Padrão | Significado |
|---|---|---|
| `--population` | o único gate final, se houver um | o caminho do gate, como `Cells/Live/CD45` |
| `--mode` | `equal_per_sample` | `equal_per_sample`, `equal_per_group` ou `ceiling` |
| `--n` | os dados decidem | por amostra, ou por grupo em `equal_per_group`; obrigatório em `ceiling` |
| `--min-events` | 1 | amostras com menos ficam de fora |
| `--group-by` | `Group` | coluna(s) que definem um grupo, separadas por vírgula |
| `--seed` | 0 | do sorteio |
| `--dry-run` | | só imprime o plano |
| `--normalize-by` | | a coluna que diz o lote de cada amostra |
| `--reference` | | as amostras que deveriam ser iguais em todo lote, como `Group=Control` |
| `--normalize-clusters` | 1 | ajusta dentro de tantos clusters do FlowSOM; 1 alinha cada marcador inteiro |
| `--normalize-to` | a média | um lote para o qual os outros são trazidos, deixado como adquirido |
| `--out` | `concat` | grava `NOME.npz` e `NOME.fcs` (e `NOME_qc.csv` com `--qc`) |

Ele imprime o plano — cada amostra com o seu grupo, o que tem e o que cede — os totais por grupo
e a fração aproveitada. Veja [[concatenation]] e [[normalisation]].

## analyze: UMAP, tSNE, FlowSOM, PhenoGraph, PCA

```
milc analyze concat.npz --channels CD8,CD4,CD19 [--umap] [--tsne] [--flowsom]
     [--phenograph] [--pca] [--out analysis]
```

Recebe um `.npz` do `concat`, ou qualquer `.fcs` concatenado (inclusive o do FlowJo;
`--sample-column`, padrão `SampleID`, indica a coluna que diz a amostra). Sem `--channels` ele
lista os canais e para. Os canais podem ser dados pelo detector ou pelo marcador.

| Opção | Padrão | Significado |
|---|---|---|
| `--transform` | `biex` | `biex` (o que o FlowJo passa aos plugins), `arcsinh` ou `logicle` |
| `--cofactor` | 150 | para `arcsinh` |
| `--neighbors`, `--min-dist` | 15, 0,5 | UMAP |
| `--perplexity` | 30 | tSNE |
| `--clusters`, `--grid` | 5, 10 | metaclusters do FlowSOM e lado da grade |
| `--k`, `--resolution` | 30, 1,0 | vizinhos do PhenoGraph; acima de 1, clusters mais numerosos e menores |
| `--method` | `leiden` | ou `louvain` |
| `--min-cluster-size` | 10 | clusters do PhenoGraph menores se juntam ao vizinho mais ligado |
| `--seed`, `--group-by` | 0, `Group` | |

Ele grava `NOME.npz` (o conjunto com os resultados) e `NOME_events.csv` (cada resultado por
evento), e para cada agrupamento as medianas por cluster (`_mfi.csv`), as frequências por amostra
(`_freq.csv`) e, havendo grupos, a média e o DP delas (`_freq_mean.csv`, `_freq_sd.csv`); a PCA
grava os seus loadings. Veja [[embeddings]] e [[clustering-methods]].

## stats: a tabela de estatísticas

```
milc stats DADOS [--wsp ARQUIVO] [--population "A;B"] [--channels CD8,CD4]
     [--stats Median,CV] [--out statistics.csv] [--prism ARQUIVO.pzfx]
```

Uma linha por população de cada amostra: contagens, porcentagens, e as estatísticas escolhidas —
`Median`, `Mean`, `Geometric mean`, `SD`, `CV`, `Robust SD`, `Robust CV`, `MAD` — de cada canal.
`--prism` grava também um arquivo do Prism com uma coluna (`--prism-value`, padrão `% of
parent`), agrupada por `--group-by`, com uma tabela para cada valor de `--within`. Veja
[[statistics-tab]] e [[tables-and-prism]].

## export: os eventos de uma população

```
milc export DADOS [--population CAMINHO] [--format fcs|csv] [--raw] [--out exported]
```

Um arquivo por amostra. `--raw` grava os valores como adquiridos, sem compensar; um arquivo FCS
leva então a matriz.

## gates: os gates para outro programa

```
milc gates DADOS [--wsp ARQUIVO] --out gates.wsp      # ou gates.xml
```

Um workspace do FlowJo se o nome terminar em `.wsp`, senão Gating-ML 2.0 (`--format
wsp|gatingml` passa por cima do nome). `--no-counts` evita abrir os arquivos para gravar a
contagem de cada população. Sem `--wsp`, os gates vêm do workspace da pasta ou dos salvos ao lado
dos arquivos. Veja [[workspaces-and-gating-ml]].

## qc: a qualidade de cada aquisição

```
milc qc DADOS [--out qc.csv] [--plots qc.pdf] [--plot-format png|pdf] [--plot-signals 3]
```

Para cada arquivo, o que seria removido e o que achou; nada é removido. `--plots ARQUIVO.pdf`
grava um PDF com uma página por amostra; qualquer outro caminho é uma pasta que recebe um arquivo
por amostra. Veja [[quality-control]].

Os limites, no `qc` e onde quer que `--qc` seja aceito: `--mad` (6), `--consecutive` (5),
`--events-per-bin` (auto), `--drift` (0,25), `--shift` (0,25), e `--qc-channels` (no `qc`
também `--channels`) para os parâmetros a julgar. `--qc` em `counts`, `concat`, `stats`,
`export` e `gates` remove os trechos ruins só naquela execução: Clean vira a primeira população,
os caminhos dos gates são escritos como antes, e o projeto ao lado dos arquivos não muda. Uma
pasta em que o controle de qualidade foi aplicado na janela é usada como está.

## autocomp: uma matriz a partir de controles de uma cor

```
milc autocomp CONTROLES… [--method autospill|classic] [--gate] [--out spillover.csv]
```

Os arquivos de controle, ou uma pasta deles: um por fluorocromo e, opcionalmente, um não marcado.
`--assign ARQUIVO=DETECTOR` (pode repetir) define o detector de um controle, `ARQUIVO=unstained`
marca o não marcado; `--detectors`, `--universal-negative` e `--max-events` (100.000; 0 para
todos) refinam. Veja [[autocompensation]].

## selftest e gui

`milc selftest` roda o caminho inteiro em dados inventados e diz se esta instalação está sã
(`--no-window` pula a janela, `--log ARQUIVO` grava num arquivo e nunca levanta erro); veja
[[installation]]. `milc gui [PASTA]` abre a janela, com uma pasta se for dada.
