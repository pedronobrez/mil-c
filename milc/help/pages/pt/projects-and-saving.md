---
title: Projetos e salvamento
---

O MIL-C não tem um comando **Save project**. O que você faz com uma pasta é gravado, à medida que
você faz, num arquivo ao lado dos `.fcs`; abrir a pasta de novo traz tudo de volta. Os eventos
concatenados, que são a parte cara de refazer, são salvos à mão.

## O arquivo do projeto: MIL-C.json

O `MIL-C.json` fica na pasta dos arquivos `.fcs`. Ele é gravado sempre que você muda algo —
digita na tabela da aba Samples, marca uma amostra, desenha um gate, escolhe uma escala — e nunca
só porque a pasta foi aberta: uma pasta apenas olhada fica como estava. Ele guarda:

| O quê | Definido em |
|---|---|
| as colunas de anotação, seus valores para cada arquivo, e quais colunas formam os grupos | [[samples-tab]] |
| quais arquivos estão desmarcados | [[samples-tab]] |
| a população escolhida para concatenar | [[samples-tab]] |
| a ordem dos valores de cada coluna | [[sample-order]] |
| os ajustes da normalização de lotes | [[normalisation]] |
| os ajustes do controle de qualidade, e quais amostras foram verificadas | [[quality-control]] |
| os gates de cada amostra cujos gates foram desenhados ou mudados aqui | [[drawing-gates]] |
| a escala escolhida para cada parâmetro | [[scales-and-transforms]] |
| uma matriz de compensação editada, comum ou por amostra | [[compensation]] |
| o unmixing espectral | [[spectral-unmixing]] |
| se os gates ainda precisam de uma olhada depois que os valores mudaram | [[gates-tab]] |
| as pranchas de painéis | [[the-plate]] |

É JSON simples: legível, e pequeno o bastante para ir num controle de versão ou junto com os
dados. Os arquivos são citados nele pelo nome, então a pasta pode ser movida ou copiada inteira.

As pranchas guardam os seus arrays num segundo arquivo, `MIL-C plates.npz`, ao lado do primeiro.
Mantenha os dois juntos. Uma pasta analisada antes de existir o arquivo de projeto pode ter um
`MIL-C gates.json`; ele é lido quando não há `MIL-C.json`, e a próxima mudança grava um
`MIL-C.json`.

O que não está no projeto: o workspace do FlowJo (lido, nunca alterado), os eventos concatenados
e o que foi calculado sobre eles (salve-os como abaixo), e os ajustes da seção Downsampling —
**Balance**, **N**, **Minimum** e **Seed** começam nos valores padrão numa sessão nova. A semente
usada de fato fica registrada em cada concatenação salva.

A pasta precisa permitir gravação. Uma pasta num compartilhamento só de leitura abre, mas as
mudanças não podem ser guardadas.

## Os eventos concatenados

**Concatenate** na aba Samples monta um novo conjunto de eventos na memória. Para guardá-lo, com
os embeddings e agrupamentos rodados sobre ele:

- `File ▸ Save concatenated file…` (`Ctrl+S`, `⌘S` no macOS). O diálogo propõe `concat.npz` na
  pasta e oferece dois formatos:

| Formato | Guarda | Para |
|---|---|---|
| conjunto do MIL-C (`.npz`) | os eventos, de que amostra e de que evento do arquivo veio cada um, cada resultado (UMAP, tSNE, FlowSOM, PhenoGraph, PCA), como cada um foi rodado, e como o conjunto foi montado | reabrir no MIL-C |
| FCS (`.fcs`) | os eventos, os canais compensados chamados `Comp-…`, e uma última coluna `SampleID` (1, 2, 3… na ordem das amostras) | o FlowJo ou qualquer outro programa |

- `File ▸ Open concatenated file…` lê qualquer um dos dois de volta. Um `.npz` abre com os seus
  resultados, então os mapas e clusters estão lá sem rodar nada de novo.

Sem nada concatenado, salvar diz `Nothing concatenated yet` na barra de status.

O `.npz` registra como os eventos foram tirados: a população, o modo de balanceamento, o que cada
amostra tinha e cedeu, quais ficaram de fora e por quê, a amostra que limitou, a semente e — onde
foram aplicados — os registros do controle de qualidade e da normalização de lotes. O arquivo FCS
leva os nomes das amostras (keyword `MILC SAMPLES`) e, onde houve controle de qualidade, quanto
ele tirou de cada amostra (keyword `MILC QC`).

## Outros arquivos que o MIL-C grava

| Comando | Grava |
|---|---|
| `File ▸ Export tables…` | `cluster_expression.csv` e `cluster_frequency.csv`, depois do FlowSOM ([[cluster-statistics]]) |
| `File ▸ Export gates…` | um workspace do FlowJo ou Gating-ML 2.0 ([[workspaces-and-gating-ml]]) |
| **Export…** na aba Gates | os eventos de uma população, um arquivo FCS ou CSV por amostra ([[gates-tab]]) |
| **Figure…** em qualquer gráfico | SVG, PDF, PNG ou TIFF ([[figures]]) |

Os arquivos CSV são gravados com marca de ordem de bytes (BOM), para que o Excel abra
corretamente nomes com acento.

## Quando os eventos salvos ficam desatualizados

Os eventos concatenados são uma cópia de valores tirados num momento. Se os valores por baixo
deles mudarem depois — outra matriz de compensação, outro unmixing, normalização de lotes
aplicada às próprias amostras (ou desligada de novo) — aparece uma linha acima das abas:

> The concatenated events were taken before … : Explore and Clusters still show the old values.
> Concatenate again in Samples to bring them up to date.

Ou seja: Explore e Clusters ainda mostram os valores antigos. Nada é recalculado às suas costas;
a linha fica até você apertar **Concatenate** de novo.

## Preferências do programa

O que é seu e não de uma pasta — a aparência (`View ▸ Appearance`), a última pasta aberta, o
idioma do manual — fica no lugar usual do sistema para preferências de aplicativos, sob o nome
MIL-C.
