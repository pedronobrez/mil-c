---
title: Workspaces do FlowJo e Gating-ML
---

O MIL-C lê um workspace do FlowJo 10 (`.wsp`) naquilo de que a análise precisa — quais arquivos, a matriz de compensação dada a cada um, e a árvore de gates de cada amostra com as contagens que o FlowJo registrou — e grava os gates de volta para outro programa, como workspace do FlowJo ou como Gating-ML 2.0, o padrão da ISAC. Gating-ML é gravado, não lido.

## Abrir um workspace

Um workspace na pasta é carregado quando a pasta é aberta (o primeiro `.wsp`, pelo nome, se houver vários). Outro pode ser escolhido em **Workspace…** na [[samples-tab]] ou em **File ▸ Open workspace…**; antes é preciso ter uma pasta aberta. A linha sob o nome da pasta diz de que workspace vêm os gates e a compensação.

As amostras do workspace são casadas com os arquivos da pasta pelo nome do arquivo: não importa onde o FlowJo achava o arquivo (muitas vezes, um caminho em outro computador). Um arquivo que o workspace lista recebe os seus gates e a sua matriz. Um arquivo que ele não lista não recebe gates — a menos que todas as amostras do workspace tenham os mesmos gates e a mesma matriz; nesse caso aparece **Gate the files the workspace does not list** na aba Samples. Marcado, esses arquivos recebem os gates e a matriz comuns, como faz o FlowJo quando se solta um arquivo no grupo. Ninguém olhou esses gates nesses arquivos: confira.

Os gates desenhados ou movidos no MIL-C ficam guardados com o projeto (veja [[projects-and-saving]]) e passam a valer no lugar dos do workspace. **Back to the workspace's gates** na [[gates-tab]] os esquece.

## O que é lido, e com que confiança

| Gate no FlowJo | Lido como | Conferido |
|---|---|---|
| Retângulo, faixa | retângulo | com as contagens do FlowJo |
| Polígono | polígono, testado na escala biexponencial ou Logicle do eixo, como faz o FlowJo | com as contagens do FlowJo |
| Quadrante | quatro populações de quadrante | pelas regras de outros leitores |
| And, Or, Not | combinação | pelas regras de outros leitores |
| Elipse | elipse, convertida da grade de 256 canais do FlowJo | pelas regras de outros leitores |
| Quadrante curly | quadrante curly | com o CytoML (abaixo) |

Retângulos e polígonos são lidos por regras conferidas com as contagens do próprio FlowJo num experimento real, 72 populações. O FlowJo conta numa grade por uma regra que não é reproduzida aqui, então igualdade exata não é esperada: uma contagem é apontada como diferente quando está a mais de 1,2% (ou 3 eventos) da do FlowJo, uma margem definida a partir do que se mediu nessas 72. Os outros tipos são lidos como dois outros leitores do formato os leem, o FlowKit (Python) e o CytoML com o cytolib (R), onde os dois concordam. Nenhum workspace à mão os tinha, então:

- uma população lida assim aparece em *itálico* na aba Gates, e a sua dica dá a contagem do FlowJo ao lado da daqui: comparar as duas é a conferência;
- o que não pôde ser lido com confiança não é lido. A barra de status diz quantas populações ficaram de fora e de que tipos; tudo o que está sob uma delas vai junto;
- nós de plug-in (UMAP, FlowSOM, exploradores de clusters) são pulados: os resultados deles são calculados de novo aqui, não importados.

Os gates em Time são convertidos dos segundos do FlowJo para as unidades do arquivo com o `$TIMESTEP` do arquivo.

### O quadrante curly

No quadrante curly do FlowJo, os braços de cima e da direita se curvam para longe do centro, acompanhando o espalhamento de uma população negativa compensada. O workspace guarda só o centro e qual é o canto, não a curva. O MIL-C desenha a curva como o CytoML: em 255 canais de cada eixo, o braço da direita sobe 0,001·d², sendo d os canais a partir do centro, e o braço de cima se inclina para a direita na mesma medida; os outros dois braços são retos. Contra o CytoML, em 40.000 eventos simulados, os cantos têm os mesmos eventos a menos de até 65, que vêm da tabela da biexponencial do próprio cytolib; com aquela tabela no lugar desta, difere no máximo 1 evento. **A curva não foi comparada com o próprio FlowJo**: não se conseguiu um workspace com quadrantes curly e os seus arquivos. Outro leitor ajustou uma curva diferente às contagens do FlowJo, e as duas podem diferir em um quarto de um canto pequeno; por isso, compare a contagem do FlowJo com a daqui (`docs/comparacao.md`, seção 7).

## Exportar os gates

**File ▸ Export gates…** abre **Export gates**:

| Controle | O que faz |
|---|---|
| **Format** | **FlowJo workspace (.wsp)** ou **Gating-ML 2.0 (.xml)**. A extensão do arquivo acompanha. |
| **Only the samples ticked in Samples** | Desmarcado, todas as amostras com gates. |
| **Write each population's count** | Abre todos os arquivos para gravar a contagem daqui de cada população, para que o outro programa tenha uma contagem para comparar com a sua. |
| **Into** | O arquivo; **Choose…** para procurar. |

Ao terminar, **Gates exported** diz o que foi gravado e lista o que não foi, e por quê.

### Como workspace do FlowJo

Um workspace com todas as amostras: o arquivo de cada uma, a sua matriz de compensação, a escala de cada parâmetro e os seus gates, declarando-se do FlowJo 10.8.1.

- Um retângulo é gravado como é; um quadrante, como o FlowJo grava: um retângulo aberto de dois lados para cada canto.
- Um polígono é gravado nas unidades do arquivo, na escala em que foi desenhado; se foi desenhado numa escala em que não se sabe que o FlowJo teste (arcsinh, logarítmica), ganha vértices ao longo dos lados até virarem a mesma linha em qualquer escala.
- Uma elipse é gravada nos 256 canais do FlowJo, pelos focos e por quatro pontos da borda, em eixos lineares e biexponenciais; nos outros, vira polígono.
- Um quadrante curly é gravado como o FlowJo grava, centro e canto; um quadrante com braços inclinados vira o polígono de cada canto.
- Uma combinação vira um AndNode, OrNode ou NotNode; Time é gravado em segundos.

Ele foi lido de volta pelo MIL-C e pelo FlowKit, evento a evento. Também foi aberto à mão no FlowJo antes da versão 1.0. Quando o workspace importar, compare as contagens lá: um gate para o qual o FlowJo não tem nome fica de fora (veja o relatório mostrado depois da exportação).

### Como Gating-ML 2.0

O Gating-ML não tem amostras, então há um arquivo quando todas as amostras têm os mesmos gates e a mesma matriz, e senão um arquivo por amostra, com o nome dela. A matriz é um `spectrumMatrix`; retângulos não têm transformação; polígonos são gravados na escala em que foram desenhados (`logicle`, `fasinh`, `flog`, ou nenhuma para linear). A biexponencial do FlowJo não tem forma em Gating-ML: é gravada como a Logicle mais próxima, com vértices acrescentados onde um lado se curvaria, para que o polígono contenha os mesmos eventos. Uma elipse é um `EllipsoidGate`, exceto num eixo biexponencial, onde vira polígono. Cada gate leva o seu nome, caminho e contagem num `custom_info`. O arquivo é válido pelo esquema do padrão, e o FlowKit o lê e separa os mesmos eventos. O flowCore e outros programas também o leem.

### O que não pode ser exportado

Um cluster transformado em população (veja [[clusters-tab]]) é o conjunto dos eventos mais próximos de alguns nós de um mapa que nenhum outro programa tem: fica de fora, com tudo o que está abaixo dele. A população criada pelo [[quality-control]] é gravada como faixas de Time e a união delas. Uma população que não está em nenhuma de várias é gravada como a negação da união delas.

## Na linha de comando

`milc gates DADOS --out ARQUIVO` faz o mesmo sem a janela: `ARQUIVO.wsp` para um workspace, qualquer outro nome para Gating-ML. `--wsp` tira os gates de um workspace dado, `--format` passa por cima da extensão, `--no-counts` não abre os arquivos. Veja [[command-line]].
