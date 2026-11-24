---
title: Tabelas, exportação e Prism
---

O MIL-C entrega os seus números de quatro maneiras: a tabela de estatísticas em CSV ou na área de transferência, as tabelas dos clusters em CSV, qualquer coluna como arquivo do GraphPad Prism, e os eventos de uma população como arquivos FCS ou CSV. Esta página diz o que cada um contém.

Todo CSV que o programa escreve é UTF-8 com marca de ordem de bytes. Sem a marca, o Excel no Windows lê um CSV na página de código do sistema e `CD8 · Median` chega embaralhado; com ela, Excel, Prism, R e pandas leem corretamente. Os números vão com seis algarismos significativos; um valor ausente é uma célula vazia.

## A tabela de estatísticas

Na [[statistics-tab|aba Statistics]]:

- **Export CSV…** grava a tabela como está, todas as colunas, em `statistics.csv` ou no nome que você escolher. A pasta é lembrada.
- **Copy** põe a mesma tabela na área de transferência, separada por tabulações e com a linha de cabeçalho, para colar no Excel ou no Prism.

A tabela é longa: uma linha por amostra e população, com as colunas das amostras ao lado do nome. É o formato que R, Python e as tabelas dinâmicas do Excel querem. O Prism quer outro formato, que **Prism…** monta.

A [[command-line|linha de comando]] dá o mesmo arquivo com `milc stats`, e um arquivo do Prism com `--prism`.

## Arquivos do Prism

O Prism não aceita uma tabela longa. Ele aceita tipos próprios:

- **Tabela de colunas** — uma coluna por grupo, cada amostra uma linha na coluna do seu grupo. É o tipo em que o Prism faz um teste t ou uma ANOVA de um fator.
- **Tabela agrupada** — uma linha por população ou cluster, uma coluna por grupo e, dentro de cada coluna, uma subcoluna por amostra (as réplicas).

Os grupos raramente têm o mesmo número de amostras. Uma tabela agrupada tem um só número de subcolunas para todas as colunas, então os grupos menores são completados com células em branco; um valor ausente também fica em branco. Cada tabela abre com casas decimais suficientes para três algarismos significativos do seu menor valor, de modo que um cluster com 0,012% de uma amostra não apareça como 0,00; o arquivo guarda todos os dígitos de qualquer forma.

### Da aba Statistics

**Prism…** abre um diálogo chamado "A table for Prism":

| Campo | O que faz |
|---|---|
| **Values** | A coluna da tabela que vai para o arquivo. Começa em % of parent |
| **Columns** | As colunas das amostras que formam as colunas do Prism. Com várias marcadas, cada combinação é uma coluna |
| **A table for each** | Opcional: uma tabela para cada valor desta coluna (o dia, por exemplo), com o nome dele |

**Save…** pede o arquivo (`statistics.pzfx` por padrão). Com várias populações na tabela, cada tabela é agrupada: uma linha por população, uma coluna por grupo, uma subcoluna por amostra; uma amostra mantém a sua subcoluna ao longo das linhas, em branco onde não tem valor. Com uma só população, cada tabela é uma tabela de colunas com o nome do valor e da população, sem brancos. A tabela de estatísticas precisa ter sido calculada antes. O Prism faz os testes ele mesmo, então o diálogo não tem teste.

### Da aba Clusters

**Prism…** na [[clusters-tab|aba Clusters]] grava a fração de cada cluster em cada amostra como tabelas agrupadas chamadas "Cluster frequency (% of sample)": clusters nas linhas, grupos nas colunas, uma subcoluna por amostra. Com uma coluna em **within each** há uma tabela por valor, com o nome dele, e todas as tabelas têm as mesmas colunas, mesmo onde um grupo não tem amostra, para que fiquem alinhadas no Prism. O nome padrão é `cluster_frequency.pzfx`.

### O que foi conferido

O `.pzfx` segue os arquivos do pacote R `pzfx` (os mesmos elementos na mesma ordem), exceto que uma tabela agrupada é marcada como tabela de dois fatores, como o próprio Prism grava. Os arquivos fazem a ida e volta pelo leitor do próprio MIL-C. Foram abertos à mão no GraphPad Prism antes da versão 1.0.

A ordem dos grupos, aqui e em todas as outras tabelas, é a definida em [[sample-order]].

## As tabelas dos clusters

**File ▸ Export tables…** pede uma pasta e grava dois arquivos sobre os clusters do FlowSOM:

| Arquivo | Linhas | Colunas | Valores |
|---|---|---|---|
| `cluster_expression.csv` | Clusters, pelo nome | Os marcadores em que o FlowSOM rodou | A mediana de cada marcador em cada cluster, sobre todos os eventos concatenados, na escala em que a clusterização rodou (0 a 1023 na biexponencial) |
| `cluster_frequency.csv` | Clusters, pelo nome | Amostras | Porcentagem dos eventos concatenados da amostra no cluster; cada coluna soma 100 |

Os clusters são os que a aba Clusters mostra: juntados e renomeados. Essas tabelas são do FlowSOM; quando só o PhenoGraph foi rodado, o menu responde "Run FlowSOM first". As mesmas duas tabelas saem do `milc analyze` como `*_mfi.csv` e `*_freq.csv`.

## Os eventos de uma população

**Export…** na [[gates-tab|aba Gates]] grava os eventos da população mostrada, um arquivo por amostra. O seu diálogo, "Export *população*":

| Campo | O que faz |
|---|---|
| **Format** | **FCS 3.1** ou **CSV** |
| **Every sample ticked in Samples** | Marcado: todas as amostras marcadas na aba Samples. Desmarcado: só as amostras selecionadas em Gates |
| **Compensated values** | Marcado (o padrão): a fluorescência como este programa trabalha com ela, em colunas chamadas `Comp-…`, para que nada adiante a compense de novo. Desmarcado: os eventos como adquiridos; um arquivo FCS leva então a matriz em uso como `$SPILLOVER`, para quem o ler depois aplicar |
| **Into** | A pasta: por padrão uma pasta chamada `exported` dentro da pasta dos dados; **Choose…** para escolher outra |

Cada arquivo se chama `amostra_população.fcs` (ou `.csv`), por exemplo `sample_01_CD4.fcs`. Um CSV tem uma linha de cabeçalho com os nomes dos marcadores (o do canal onde não há marcador) e uma linha por evento. Um arquivo FCS mantém `$CYT`, `$DATE`, `$BTIM`, `$ETIM` e `TUBE NAME` do original, e acrescenta `MILC POPULATION` e `MILC SOURCE` (o arquivo de onde veio). Uma amostra sem essa população é pulada.

Os gates em si, e não os seus eventos, saem por **File ▸ Export gates…**: ver [[workspaces-and-gating-ml]].
