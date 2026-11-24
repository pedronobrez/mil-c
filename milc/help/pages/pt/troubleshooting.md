---
title: Solução de problemas
---

O que dá errado com mais frequência, o que o MIL-C diz quando acontece, e o que fazer. As
mensagens são citadas como o programa as escreve, em inglês; a maioria aparece na barra de status,
embaixo da janela, ou numa linha na cor de aviso dentro da aba ou do diálogo em questão. Um erro
num trabalho feito em segundo plano (uma concatenação, um UMAP) aparece na barra de status como o
tipo do erro seguido da mensagem, por exemplo `ValueError: …`.

## A primeira abertura

**O macOS diz que o aplicativo não pode ser aberto, ou que é de um desenvolvedor não
identificado.** Os instaladores não são assinados. Clique com o botão direito (ou Control-clique)
em **MIL-C** em Aplicativos, escolha **Abrir**, e **Abrir** de novo no diálogo. Se o macOS só
oferecer **Mover para o Lixo**, feche, abra **Ajustes do Sistema ▸ Privacidade e Segurança**,
desça até a linha sobre o MIL-C e escolha **Abrir Mesmo Assim**. Isso é pedido uma vez. Num
terminal, `xattr -dr com.apple.quarantine /Applications/MIL-C.app` faz o mesmo.

**O Windows diz que protegeu o computador (SmartScreen).** Escolha **Mais informações** e depois
**Executar assim mesmo**. Onde nenhum instalador é permitido, use o `.zip`: descompacte e rode
`MIL-C.exe`.

**Não abre num Mac Intel.** A versão para macOS é só para Apple Silicon; rode a partir do
código-fonte ([[installation]]).

**Não abre num Linux mais antigo.** O pacote precisa de uma biblioteca C tão recente quanto a do
Ubuntu 22.04. Instale pelo `.deb`, que traz as bibliotecas do sistema de que o Qt precisa, ou rode
a partir do código-fonte.

**O primeiro UMAP ou FlowSOM demora muito.** O código numérico é compilado da primeira vez e
guardado: medido no pacote do macOS, 16 s para o primeiro UMAP, 5 s para cada um depois.

**Esta instalação está sã?** Rode `MIL-C selftest` ([[command-line]]). No Windows acrescente
`--log selftest.log` e leia o arquivo: um programa de janela não tem console.

## Abrindo arquivos

| Mensagem | O que fazer |
|---|---|
| `No .fcs files in` a pasta | os arquivos precisam estar diretamente na pasta escolhida, não em subpastas |
| Note `cannot be read: not an FCS file` | o arquivo está danificado ou não é FCS; ele fica de fora e o resto abre |
| Note `cannot be read: $MODE …: only list mode is supported` | arquivos FCS em modo histograma não são lidos |
| Note `cannot be read: $DATATYPE … is not supported` | um formato numérico incomum; exporte o arquivo de novo do programa de aquisição |
| Note `different channels from the others` | os canais do arquivo diferem dos da maioria — muitas vezes é a exportação concatenada de outro programa. Tire-o da pasta |
| `Open a folder first` | um workspace é lido para uma pasta: abra a pasta, depois **Workspace…** |
| Note `no gates in the workspace` | o workspace não cita este arquivo. Marque **Gate the files the workspace does not list**, ou escolha **All events** |
| Status `N population(s) not read (…)` depois de um workspace | esses gates são de um tipo que o MIL-C não lê; o que está abaixo deles também falta |

Populações em itálico na aba Gates foram lidas do workspace por regras não conferidas com o
FlowJo; a contagem do próprio FlowJo está no tooltip para comparar
([[workspaces-and-gating-ml]]).

## Concatenando

**Concatenate está apagado.** O plano não tira nenhum evento: nenhuma amostra marcada, a população
não existe em nenhuma amostra marcada, ou as contagens ainda estão rodando (`counting…` na coluna
Note). Leia a coluna Note e o resumo abaixo dos ajustes ([[samples-tab]]).

**Uma amostra impõe um número muito baixo a todas.** A Note dela diz `sets the limit`. Suba
**Minimum** acima do tamanho dela para deixá-la de fora, desmarque-a, ou balanceie por grupo
([[concatenation]]).

**`… look like the same animal acquired twice: untick one, or it counts double.`** Dois arquivos
marcados têm nomes que diferem só por um sinal no fim, o jeito usual de nomear um tubo
readquirido. Desmarque o que você não quer.

**Uma linha acima das abas diz que os eventos concatenados foram tirados antes de algo mudar.** A
compensação, o unmixing ou a normalização de lotes mudou depois que você concatenou. Explore e
Clusters ainda mostram os valores antigos até você apertar **Concatenate** de novo.

## Normalização de lotes

| Mensagem | O que fazer |
|---|---|
| `Needs samples with Group = Control in at least two batches; there are such samples in 1.` | o valor de referência falta em todos os lotes menos um, ou essas amostras estão desmarcadas ou sem a população |
| `Batch 3, which the others are to be brought onto, has no sample with Group = Control.` | escolha outro **Bring to**, ou a média dos lotes |
| `… batch-and-cluster pair(s) had too few reference events and were left as they were` | baixe **Within**, ou use mais eventos de referência: são precisos 50 por lote em cada cluster |
| `Not normalised, their batch having no reference: …` | o lote dessas amostras não tem amostra de referência |

Veja [[normalisation]].

## Controle de qualidade

`Cannot be checked:` seguido do motivo — `the file has no Time parameter` (sem parâmetro Time),
`… events are too few (it takes …)` (eventos de menos), `Time is not a number for some events`,
`Time is the same for every event`, `Time runs backwards: several acquisitions in one file, or a
clock that was reset`. Um arquivo assim fica inteiro. `That is more than a fifth of the file: look
at the plots before applying.` é um aviso, não um erro: algo grande deu errado, ou os limites estão
apertados demais para este arquivo. Veja [[quality-control]].

## Gates e compensação

- `That gate needs two parameters: choose a Y axis` — num histograma só dá para desenhar uma
  faixa.
- `Every cell has to be a number` / `This matrix cannot be inverted: …` — no diálogo de
  compensação, uma célula ficou vazia, ou duas linhas são iguais ([[compensation]]).
- `These files carry no compensation matrix` — depois de **Compensation…**: os arquivos não têm a
  keyword `$SPILL` e nenhum workspace lhes dá uma matriz, então não há o que editar. Leia um
  workspace que tenha uma ([[compensation]]); `milc autocomp` calcula uma matriz a partir de
  controles de uma cor ([[autocompensation]]).
- Uma barra na aba Gates diz que os gates foram desenhados sobre valores que mudaram desde então.
  Olhe-os e aperte **Reviewed** ([[gates-tab]]).
- `No gate proposed: …` — o gate automático não achou o que propor; o motivo vem em seguida
  ([[auto-gates]]).

## Salvando e exportando

| Mensagem | Significado |
|---|---|
| `Nothing concatenated yet` | salvar exige uma concatenação antes |
| `Run FlowSOM first: the tables are per cluster` | `File ▸ Export tables…` exige o FlowSOM |
| `No sample has gates to export` | desenhe gates ou leia um workspace antes |
| `Nothing on the plot yet` | **Figure…** exige um gráfico |
| `cannot write .xyz: use .svg, .pdf, .png, .tiff` | uma figura é gravada num desses formatos |

As mudanças são gravadas em `MIL-C.json` na pasta dos dados. Se a pasta for só de leitura (um
drive compartilhado, um CD), copie-a antes para um lugar com permissão de gravação
([[projects-and-saving]]).

## Se nada disso resolver

Relate o problema em `https://github.com/pedronobrez/mil-c/issues`, com a versão (`milc
--version`), o sistema, a mensagem como aparece e o que você fez antes dela. Não anexe dados que
você não pode compartilhar; um arquivo inventado que mostre o mesmo problema é melhor. Os termos
usados nas mensagens estão no [[glossary]].
