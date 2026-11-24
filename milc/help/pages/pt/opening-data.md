---
title: Abrindo os dados
---

O MIL-C trabalha sobre uma pasta: cada arquivo `.fcs` nela é uma amostra, um workspace do FlowJo
ao lado deles os gateia, e o que você faz com eles fica guardado na mesma pasta (veja
[[projects-and-saving]]). Não há etapa de importação nem cópia dos arquivos: eles são lidos onde
estão, e nunca são alterados.

## Abrindo uma pasta

- `File ▸ Open folder…` (`Ctrl+O`, `⌘O` no macOS), ou o botão **Open folder…** no alto, à
  direita, da aba Samples. O diálogo começa na pasta aberta da última vez.
- Pela linha de comando, `milc gui caminho/da/pasta` abre a janela já com a pasta aberta.

Todo arquivo terminado em `.fcs` diretamente dentro da pasta entra na lista; subpastas não são
vasculhadas. Os arquivos são ordenados como uma pessoa conta — `3Ctrl1` antes de `14Ctrl1`,
`sample_2` antes de `sample_10` — e de início só os cabeçalhos são lidos: o número de eventos,
os canais e os marcadores. Os eventos em si são lidos quando a amostra é necessária pela
primeira vez, para uma contagem, um gráfico ou uma concatenação. A barra de status diz quantos
arquivos foram achados, ou `No .fcs files in` a pasta.

O MIL-C lê arquivos FCS 2.0, 3.0 e 3.1 em modo lista, o formato que todo citômetro grava.

## Arquivos que ficam de fora

Dois tipos de arquivo aparecem na tabela da aba Samples mas não podem entrar. A caixa deles fica
desmarcada e apagada, e a coluna **Note** diz por quê:

| Note | Por quê |
|---|---|
| `cannot be read: …` | o cabeçalho ou os dados não puderam ser lidos; o resto da mensagem é o motivo, como `not an FCS file` ou `$DATATYPE … is not supported` |
| `different channels from the others` | os arquivos de um experimento têm os mesmos canais; um que não tem — a exportação concatenada de outro programa, por exemplo — fica de lado. Vale o conjunto de canais que a maioria dos arquivos tem |

Um arquivo quebrado nunca impede a pasta de abrir.

## O primeiro palpite de grupos

Uma pasta nova começa com uma coluna de anotação, **Group**, preenchida por um palpite tirado dos
nomes dos arquivos: o trecho de letras que assume poucos valores diferentes entre os arquivos.
`Control_1`, `Control_2`, `Treated_1`, `Treated_2` dão Control e Treated. Onde nenhum trecho de
letras serve, a coluna começa vazia. O palpite é só um ponto de partida: digite por cima,
acrescente colunas ou leia-as dos nomes com um padrão, como [[samples-tab]] descreve.

## O workspace do FlowJo

Se a pasta tiver um arquivo `.wsp`, o primeiro em ordem alfabética é lido junto com a pasta.
Para usar outro, ou acrescentar um depois, use **Workspace…** na aba Samples ou `File ▸ Open
workspace…`; é preciso ter uma pasta aberta antes, senão a barra de status diz `Open a folder
first`.

Do workspace o MIL-C toma, para cada arquivo que ele cita (casado pelo nome do arquivo):

- a árvore de gates, guardando a contagem do FlowJo de cada população para comparar;
- a matriz de compensação que aquela amostra recebeu.

Os eixos são desenhados nas escalas do próprio MIL-C — a biexponencial do FlowJo para
fluorescência, os valores como gravados para dispersão e tempo — até você escolher outras
([[scales-and-transforms]]).

Retângulos, polígonos, faixas, quadrantes, combinações booleanas, elipses e os quadrantes curvos
("curly") do FlowJo são lidos. Os quatro últimos são lidos por regras tiradas do FlowKit e do
CytoML, não conferidas com o FlowJo em arquivos reais; as populações lidas assim aparecem em
itálico na aba Gates, com a contagem do próprio FlowJo no tooltip. Depois de ler um workspace, a
barra de status relata o que vale saber: populações de um tipo que não pôde ser lido, e quantas
foram lidas por regras não conferidas. [[workspaces-and-gating-ml]] entra nos detalhes.

Abaixo do nome da pasta, a aba Samples diz o que está em uso: `Gates and compensation from` o
workspace, ou, sem ele, `No workspace: every event of every file, compensated with the file's own
matrix` (todos os eventos de cada arquivo, compensados com a matriz do próprio arquivo).

## Arquivos que o workspace não cita

Muitas vezes o workspace cita só parte dos arquivos da pasta — um tubo adquirido depois que a
análise foi montada, por exemplo. Quando todas as amostras do workspace têm os mesmos gates e a
mesma matriz de compensação, aparece na aba Samples a caixa **Gate the files the workspace does
not list**. Marcada (o padrão), os outros arquivos recebem esses gates e essa matriz, que é o que
o FlowJo faz quando se solta um arquivo num grupo; a Note deles diz `gates from the group`.
Ninguém olhou esses gates nesses arquivos, então olhe-os em [[gates-tab]]. Desmarcada, eles
ficam sem gates e só podem entrar quando **Population** é **All events**; com uma população
gateada, a Note diz `no gates in the workspace`.

## O que é concatenado por padrão

A **Population** da aba Samples começa no gate mais profundo que toda amostra gateada tem — o fim
da cadeia mais longa, como `Cells/Single/Live/CD45`, e não um gate desenhado ao lado dela. A
escolhida da última vez é lembrada com o projeto, desde que os gates ainda a tenham.

## Abrindo um arquivo concatenado

`File ▸ Open concatenated file…` abre eventos concatenados antes, pulando a aba Samples:

- um `.npz` salvo pelo MIL-C (`File ▸ Save concatenated file…`), com tudo o que foi calculado
  sobre ele;
- um `.fcs` concatenado, do MIL-C ou do FlowJo. Canais chamados `Comp-…` ou `FJComp-…` são
  tomados como compensados, e uma coluna `SampleID`, quando existe, diz de que amostra veio cada
  evento.

A janela passa então para a aba Explore. Veja [[concatenation]] e [[embeddings]].

## Em seguida

Com a pasta aberta, a tabela de [[samples-tab]] é onde os arquivos são agrupados e escolhidos;
[[quality-control]] pode limpar cada aquisição antes de tirar eventos dela.
