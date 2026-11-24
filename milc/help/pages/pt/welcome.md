---
title: Bem-vindo ao MIL-C
---

O MIL-C é um programa de computador para citometria de fluxo. Ele recebe uma pasta de arquivos
`.fcs` e, quando houver, o workspace do FlowJo que os gateia; tira o mesmo número de células de
cada amostra; roda UMAP, tSNE, FlowSOM e um agrupamento no estilo do PhenoGraph sobre o que tirou;
compara grupos; e transforma qualquer gráfico numa figura no tamanho, na fonte e nas cores que a
revista pede.

Ele foi escrito para substituir um caminho no FlowJo — gatear, DownSample, concatenar e depois os
plugins de UMAP, FlowSOM e ClusterExplorer — que funcionava, mas deixava três coisas a desejar:
figuras difíceis de acertar, uma concatenação que favorecia em silêncio a amostra com mais
células, e nenhum registro das sementes e parâmetros por trás de um resultado. O MIL-C é código
aberto (GPL-3.0-or-later), gratuito e não precisa de R instalado. Ele lê workspaces do FlowJo e
foi conferido com os números do FlowJo, mas não contém código do FlowJo e não tem vínculo com o
FlowJo nem com a BD.

## A janela

A janela principal tem cinco abas, lidas da esquerda para a direita à medida que a análise
avança. Os nomes ficam em inglês, como na tela:

| Aba | Para que serve |
|---|---|
| **Samples** | quais arquivos entram, como se agrupam, quantos eventos cada um cede, normalização de lotes |
| **Gates** | uma amostra por vez: gráficos, a árvore de gates, compensação, escalas, qualidade, modelos |
| **Statistics** | uma linha por população de cada amostra, testes entre grupos, arquivos do Prism |
| **Explore** | UMAP, tSNE, FlowSOM, PhenoGraph e PCA sobre os eventos concatenados |
| **Clusters** | os clusters: mapa, heatmap, árvore e frequências, ligados; abundância e expressão |

## O manual, seção por seção

| Seção | Páginas |
|---|---|
| Primeiros passos | [[installation]] · [[opening-data]] · [[projects-and-saving]] |
| Amostras (Samples) | [[samples-tab]] · [[concatenation]] · [[sample-order]] · [[normalisation]] · [[quality-control]] |
| Gates | [[gates-tab]] · [[drawing-gates]] · [[auto-gates]] · [[scales-and-transforms]] · [[compensation]] · [[autocompensation]] · [[spectral-unmixing]] · [[cell-cycle]] · [[proliferation]] · [[kinetics]] |
| Estatística (Statistics) | [[statistics-tab]] · [[abundance]] · [[differential-expression]] · [[tables-and-prism]] |
| Explorar (Explore) | [[explore-tab]] · [[embeddings]] |
| Clusters | [[clusters-tab]] · [[clustering-methods]] · [[cluster-statistics]] |
| Referência | [[figures]] · [[the-plate]] · [[workspaces-and-gating-ml]] · [[command-line]] · [[keyboard-shortcuts]] · [[measured-facts]] · [[design-principles]] · [[troubleshooting]] · [[glossary]] · [[version-history]] |

Numa primeira sessão costuma-se ler [[installation]], [[opening-data]], [[samples-tab]] e
[[concatenation]], depois [[explore-tab]] e [[clusters-tab]].

## Como chegar ao manual

- **Help ▸ Manual** abre esta janela. O atalho é `F1` no Windows e no Linux e `⌘?` no macOS.
- O manual abre na página de onde você está: a aba visível ou o diálogo aberto na frente. Aperte
  `F1` no diálogo de qualidade e você cai em [[quality-control]]; na aba Samples, em
  [[samples-tab]].
- **Help ▸ Export manual as PDF…** grava todas as páginas, na ordem do sumário, num único
  documento A4. Os links entre páginas continuam sendo links dentro do PDF. O botão **Export as
  PDF…** na barra do próprio manual faz o mesmo.

## Como ler

- O sumário fica à esquerda; a página, à direita. **◀ Back** e **Forward ▶** percorrem as páginas
  que você abriu, e **Contents** volta para esta.
- Um link no texto abre outra página. No pé de cada página, **Recebe links de** lista as páginas
  que apontam para ela, para você voltar pelo caminho que um leitor teria feito.
- A caixa de busca (`Ctrl+F`, `⌘F` no macOS) procura em todas as páginas enquanto você digita.
  Cada palavra casa pelo começo — `compens` acha compensação e compensado — e uma página só
  aparece se tiver todas as palavras digitadas, a melhor primeiro; uma palavra no título da página
  pesa mais. `Return` abre o primeiro resultado.
- O seletor **English / Português** na barra troca o idioma do manual, e a escolha é lembrada. O
  programa em si continua em inglês; por isso as páginas em português mantêm em inglês os nomes de
  botões e menus.

## Convenções

- **Negrito** é um controle da tela, escrito exatamente como aparece: **Concatenate**,
  **Order…**, **Apply to every sample**.
- `Código` é algo que se digita, ou o nome de um arquivo: `MIL-C.json`, `milc concat`.
- Um caminho por menus é escrito `File ▸ Open folder…`.
- Os exemplos usam nomes inventados — arquivos `sample_01.fcs`, grupos Control e Treated.
- Onde um número foi medido (contra o FlowJo, um pacote de R, dados simulados), a página diz isso
  e onde; as medidas reunidas estão em [[measured-facts]].

## O que o MIL-C não é

- Não é o FlowJo inteiro. Os gates são retângulos, polígonos, laços, elipses, faixas, quadrantes
  e combinações booleanas; o editor de layout do FlowJo, layouts de placa de 96 poços e a edição
  de keywords não estão aqui.
- Não é um programa de aquisição. Ele lê arquivos FCS 2.0, 3.0 e 3.1 gravados por um citômetro;
  não controla um.
- Não é um pacote de estatística. Os testes comparam um valor por amostra entre grupos, pareados
  ou não; modelos mistos ficam para outro programa, e [[tables-and-prism]] diz como levar os
  números até lá.
- Não é o FlowJo ponto a ponto. Os mapas e clusters ficam próximos dos do FlowJo no que dá para
  medir, não idênticos; [[design-principles]] explica as escolhas.
