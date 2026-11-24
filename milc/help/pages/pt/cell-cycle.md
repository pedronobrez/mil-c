---
title: Ciclo celular
---

O ciclo celular é lido do conteúdo de DNA de cada célula (iodeto de propídio, DAPI, Hoechst…) numa escala linear. As células com uma dose de DNA formam um pico (G0/G1), as com duas formam um pico em cerca do dobro do sinal (G2/M), e as que estão copiando o DNA (S) ficam entre os dois. Um modelo do histograma diz quantas células há em cada fase.

## Onde fica

Na [[gates-tab]], escolha a população (células únicas, a menos que os agregados sejam modelados), selecione as amostras na lista e depois **Model…** ▸ **Cell cycle…**. O diálogo ajusta todas as amostras selecionadas e mostra uma por vez.

Os valores são os que os arquivos guardam, compensados, sem escala de exibição. O **Parameter** começa no parâmetro de X se ele for linear; senão, no primeiro parâmetro linear.

## Os ajustes

| Controle | Padrão | O que faz |
|---|---|---|
| **Parameter** | | O sinal de DNA. |
| **Model** | Watson pragmatic | **Watson pragmatic** ou **Dean-Jett-Fox**, abaixo. |
| **G2 / G1** | 1,90 a 2,10 | A média de G2 fica entre estes múltiplos da média de G1: um corante raramente é exatamente proporcional ao DNA. O mesmo número duas vezes fixa a razão. |
| **G1 near** | Found | Mais ou menos onde está o pico G1, nas unidades do arquivo, quando o achado não é ele. Com uma segunda população: onde está o G1 diploide. |
| **G2 has the CV of G1** | desligado | Desmarcado, cada pico tem o seu CV. Marcado, um CV para os dois: mais estável quando o G2 é pequeno. |
| **Second (aneuploid) population** | desligado | Um segundo G1, S e G2 (só no Dean-Jett-Fox). |
| **DNA index near** | Found | Mais ou menos onde está o G1 aneuploide, como múltiplo do G1 diploide (1,5: metade a mais de DNA). |
| **Debris** | None | **None**, **Exponential** ou **Single cut** (só no Dean-Jett-Fox). |
| **Aggregates (doublets, triplets)** | desligado | Duas ou três células grudadas (só no Dean-Jett-Fox). |

O ajuste roda de novo um instante depois de qualquer mudança.

## Os modelos

Nos dois, G0/G1 e G2/M são gaussianas. O histograma tem 256 classes e é ajustado de 4 desvios-padrão abaixo do G1 a 4 acima do G2. O número de células de cada componente é achado por mínimos quadrados não negativos; só as posições e larguras são procuradas. O ajuste é por mínimos quadrados nas contagens, sem pesos, como os do FlowJo e do ModFit.

- **Watson pragmatic** (segundo Watson, Chambers e Smith, 1987) não supõe nada sobre a forma de S. O G1 é ajustado no lado do pico que S não alcança (o esquerdo), o G2 no seu lado direito, e S é o que sobra entre as duas médias depois de tirar os dois. O S que fica sob os lados ajustados é levado em conta, e os picos são ajustados de novo quatro vezes com ele. É menos preciso, porque usa metade de cada pico, e enviesado onde o S se acumula junto a um pico.
- **Dean-Jett-Fox** (Dean e Jett, 1974; Fox, 1980) dá uma forma ao S: um polinômio de segundo grau do conteúdo de DNA entre as duas médias, nunca negativo, alargado pelo mesmo erro de medida dos picos. Os três são ajustados de uma vez ao histograma inteiro. É o mais preciso quando o S é suave, e errado quando não é.

Com o Dean-Jett-Fox, mais três componentes podem ser ajustados junto com o resto:

- **Uma segunda população**: seu próprio G1, S e G2, com o G1 no índice de DNA vezes o do primeiro (procurado a até 10% de onde foi achado ou indicado), o G2 no mesmo múltiplo do seu G1 que o do primeiro, e um CV próprio. O G1 mais baixo é tomado como diploide, a menos que **G1 near** diga outra coisa. Um segundo G1 perto de 2 não pode ser distinguido do G2 diploide: nesse caso, indique **DNA index near**.
- **Debris**: **Exponential**, uma curva descendente cuja taxa é ajustada, ou **Single cut**, núcleos cortados uma vez num ponto ao acaso, calculado a partir das células intactas achadas (segundo Bagwell et al., 1991). O ajuste passa então a ir até um quinto do G1.
- **Aggregates**: dubletos e tripletos, a distribuição das células únicas convoluída consigo mesma, como faz o modelo de agregados do ModFit. G1+G1 cai sobre o G2; G1+G2 e G2+G2, acima dele. Com isso, não é preciso fazer gate de células únicas.

Debris e agregados tiram a forma das células achadas, então o ajuste é repetido quatro vezes, cada uma com as formas que a anterior deu.

## A tabela

Uma linha por amostra, com as suas colunas da [[samples-tab]]. As colunas principais:

| Coluna | Significado |
|---|---|
| **% G1**, **% S**, **% G2** | Das células da primeira população (diploide), cada uma com o seu **SE**. |
| **G1 mean**, **G1 CV**, **G2 mean**, **G2 CV** | Posição e CV de cada pico. |
| **G2/G1** | A razão entre as duas médias, com o seu SE. |
| **% sub-G1**, **% above G2** | De todos os eventos, fora da faixa ajustada (as linhas pontilhadas no gráfico). |
| **RMS**, **RMS %** | O resíduo, em eventos por classe e como porcentagem da classe mais alta. |
| **Chi²** | Qui-quadrado reduzido: perto de 1 quando só sobra o erro de contagem. |

Com os extras: **% aneuploid** (das células das duas populações), **DNA index**, as fases e o CV do aneuploide; **% debris** (de todos os eventos modelados); **% aggregates** (de eventos: um agregado é um evento).

Os erros-padrão são os da contagem: quanto as mesmas células, adquiridas de novo, mudariam o número. Eles não dizem se o modelo está certo; o RMS e o qui-quadrado dizem. Uma amostra que não pode ser ajustada tem o motivo numa coluna **Note**.

## O gráfico e os botões

O gráfico mostra o histograma, cada componente e o modelo da amostra escolhida na tabela, com uma linha acima que resume o ajuste. **Figure…** faz dele uma figura ([[figures]]); **Compare groups…** testa uma coluna entre grupos de amostras, cada amostra uma observação ([[statistics-tab]]); **Export CSV…** grava a tabela; **Copy** a põe na área de transferência para o Excel ou o Prism ([[tables-and-prism]]).

## Até onde confiar

Em histogramas de 20.000 células simuladas com fases conhecidas: o Dean-Jett-Fox fica a menos de 0,3 ponto com CV de 3 a 5% quando o S é suave, e varia de 0,5 a 1,5 ponto de uma aquisição para outra; o Watson fica a cerca de 1 ponto em média, com 1,5 a 3 pontos de dispersão. Onde o S tem uma forma que nenhuma parábola tem, o Dean-Jett-Fox pode errar de 6 a 12 pontos, e o seu qui-quadrado (de 9 a 30 em vez de 1) acusa; o do Watson não.

Em 223 arquivos públicos de BrdU/PI, em que o BrdU dá a fase de cada célula sem modelo (`docs/comparacao.md`, seção 6): nos 39 controles, o Dean-Jett-Fox deu 6,6 pontos de S a menos que o BrdU (DP 5,8), no sentido esperado, já que o BrdU também marca células que passaram por S durante o pulso; o Watson deu 3,9 a mais, com DP de 13,5. Ajustados sem gate de singletos, os agregados achados (24,8% em média) bateram com o que um gate de singletos tira (25,6%; r = 0,90). Nada foi comparado com o ModFit ou o FlowJo, que não estavam disponíveis.
