---
title: A aba Samples
---

A aba Samples é a primeira da janela e aquela de que as outras dependem. Ela tem a tabela dos
arquivos — uma linha por arquivo `.fcs` da pasta — com tudo o que você sabe de cada amostra e
que o arquivo não registra, e o painel à direita que decide quais eventos vão para o conjunto
concatenado: de que população, de quais amostras, quantos de cada, e se os lotes são alinhados
antes. Nada é tirado dos arquivos até apertar **Concatenate**; até lá, todo número da tabela é um
plano que você pode ler e mudar.

## Acima da tabela

À esquerda, a pasta aberta e, abaixo, de onde vêm os gates e a compensação: o nome do workspace,
ou `No workspace: every event of every file, compensated with the file's own matrix`. À direita,
cinco botões numa linha que quebra quando a janela é estreita:

| Botão | O que faz |
|---|---|
| **Add column…** | uma nova coluna de anotação: Day, Tissue, Donor, Treatment… |
| **Order…** | qual valor de uma coluna vem primeiro em todo gráfico e tabela ([[sample-order]]) |
| **Fill from file names…** | lê colunas dos nomes dos arquivos com um padrão (abaixo) |
| **Workspace…** | lê um `.wsp` do FlowJo e aplica os seus gates e compensação ([[opening-data]]) |
| **Open folder…** | abre outra pasta |

## A tabela

| Coluna | Significado | Editável |
|---|---|---|
| (caixa) | marcada: a amostra entra na concatenação e nas estatísticas | sim |
| **Sample** | o nome do arquivo, sem `.fcs` | não |
| **Group**, e cada coluna que você acrescenta | a sua anotação: texto livre | sim |
| **Events** | eventos no arquivo, pelo cabeçalho | não |
| **In gate** | eventos da **Population** escolhida neste arquivo | não |
| **Takes** | eventos que o plano tira desta amostra | não |
| **Note** | por que uma amostra não cede nada, ou algo a saber sobre ela | não |

Para digitar numa coluna sua, dê duplo clique na célula, ou selecione-a e comece a digitar. Com
várias linhas selecionadas (Shift-clique, ou `Ctrl`/`⌘`-clique), o que você digita numa vai para
todas: dez amostras de um grupo são uma digitação, não dez. Os espaços nas pontas dos valores são
retirados.

Dê duplo clique no cabeçalho de uma coluna sua para renomeá-la; clique com o botão direito para
**Rename…** e **Remove**. Os nomes das colunas da própria tabela — Sample, Events, In gate,
Takes, Note — não podem ser dados a uma coluna sua, e duas colunas não podem ter o mesmo nome; a
barra de status avisa se você tentar.

A coluna **Note** pode dizer:

| Note | Significado |
|---|---|
| `sets the limit` | a amostra cujo tamanho decidiu o número de todas ([[concatenation]]) |
| `gates from the group` | o workspace não cita este arquivo; ele recebeu os gates que os citados compartilham |
| `your gates` | os gates desta amostra foram desenhados ou mudados no MIL-C |
| `counting…` | a população está sendo contada neste arquivo |
| `N events, below the minimum of M` | menos eventos que **Minimum**: fica de fora |
| `no gates in the workspace` | o arquivo não tem gates, então não tem eventos de uma população gateada |
| `cannot be read: …`, `different channels from the others` | o arquivo não pode entrar ([[opening-data]]) |

## Preencher a partir dos nomes

**Fill from file names…** abre um diálogo com uma expressão regular proposta a partir do primeiro
nome de arquivo, na qual cada grupo nomeado `(?P<Nome>…)` vira uma coluna chamada Nome: um para
cada trecho de dígitos e cada trecho de letras, o primeiro trecho de letras chamado Group. Para um
arquivo `7Treated3` (dia 7, Treated, animal 3) a proposta é

```
(?P<Number1>\d+)(?P<Group>[A-Za-z]+)(?P<Number2>\d+)
```

Renomeie os grupos para o que eles são — `(?P<Day>\d+)(?P<Group>[A-Za-z]+)(?P<Animal>\d+)` — e
apague os que você não quer como coluna. `\d+` casa um número, `[A-Za-z]+` um trecho de letras.
A prévia mostra o que o padrão acha nos oito primeiros nomes; abaixo dela, quantos de todos os
nomes ele casa. Um arquivo que o padrão não casa fica com o que tinha. O **OK** fica ativo quando
o padrão tem pelo menos um grupo nomeado e casa pelo menos um nome. Colunas que ainda não existem
são criadas.

## O painel à direita

Os ajustes rolam quando a janela é baixa; **Concatenate** fica sempre visível embaixo deles.

**Samples**

- **Population**: quais eventos de cada amostra são concatenados — **All events**, ou um gate da
  árvore, recuado conforme a profundidade. Ele precisa existir em toda amostra que entra; a lista
  traz os gates que todas as amostras gateadas têm. Escolher um conta-o em cada arquivo (a barra
  de status mostra `Counting …`) e preenche **In gate**.
- **Include**: quantas amostras estão marcadas, como `8 of 12 samples`. O menu marca **Every
  sample** ou **None**, e tem um submenu por coluna (com entre 2 e 40 valores) para marcar ou
  desmarcar de uma vez todas as amostras de um valor — para concatenar só dois grupos, desmarque
  os outros.
- **Gate the files the workspace does not list**: aparece só quando se aplica; veja
  [[opening-data]].
- **Group by**: qual coluna, ou quais colunas juntas, formam um grupo — para balancear o mesmo
  total por grupo e para toda comparação depois. Marque duas e cada combinação é um grupo:
  `Treated · 7`. Sem nenhuma marcada o botão diz `Nothing: one group`.

**Downsampling** — **Balance**, **N**, **Minimum** e **Seed**, descritos em [[concatenation]].

**Batch normalisation** — **Batches**, **Reference**, **Within**, **Bring to** e **Also in Gates
and Statistics**, descritos em [[normalisation]].

## O resumo e os avisos

Abaixo dos ajustes, o plano num parágrafo: o total, o número de amostras e a fração do que elas
têm, por exemplo **72,000 events** from 6 samples — 31% of what they have; depois o total de cada
grupo e `Limited by` a amostra ou grupo que definiu o número. Na cor de aviso, abaixo:

- as notas do próprio plano (uma amostra abaixo de **N**, amostras de um grupo cedendo números
  diferentes);
- arquivos sem gates que ficam fora de uma população gateada;
- dois arquivos marcados cujos nomes diferem só por um sinal no fim — `sample_04` e
  `sample_04'` — que é como se costuma nomear a reaquisição de um tubo: `look like the same
  animal acquired twice: untick one, or it counts double` (parecem o mesmo animal adquirido duas
  vezes: desmarque um, senão ele conta em dobro).

## Concatenate

**Concatenate** executa o plano, monta o conjunto concatenado e passa para a aba Explore. Ele fica
desativado enquanto o plano não tira nada ou uma concatenação está em andamento; a barra de
status mostra `Reading amostra (3/12)` durante o processo, e `Concatenated N events from M
samples` no fim. Veja [[concatenation]] para o que é montado, e [[explore-tab]] para o que vem
depois.

Toda mudança feita nesta aba é gravada na hora no arquivo de projeto da pasta
([[projects-and-saving]]), menos os quatro ajustes de Downsampling.
