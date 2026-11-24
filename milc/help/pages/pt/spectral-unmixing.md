---
title: Unmixing espectral
---

Um citômetro espectral (Cytek Aurora e Northern Lights, Sony ID7000, BD FACSymphony A5 SE) tem mais detectores que fluorocromos. Cada fluorocromo tem um espectro: quanto dele cada detector vê, escalado para que o detector mais alto valha 1. Os valores de um evento nos detectores são tomados como uma soma dos espectros, e o unmixing acha quanto de cada fluorocromo compõe essa soma. Ele toma o lugar da [[compensation]]: com o unmixing ligado, nenhuma matriz de compensação é usada.

## Onde fica

**Unmixing…** na [[gates-tab]], ao lado de **Compensation…**. O diálogo **Spectral unmixing** tem os controles de referência à esquerda e três abas à direita.

## Os controles de referência

A lista traz todos os arquivos da pasta; **Add files…** acrescenta controles guardados em outro lugar. Os arquivos cujo nome contém `(Beads)` ou `(Cells)`, como o SpectroFlo os nomeia, começam marcados, e o fluorocromo é deduzido do nome (`A3 BUV395 (Beads)_Bead controls_1.fcs` → `BUV395`).

| Coluna | Significado |
|---|---|
| **Use** | Marque cada controle de uma cor a usar. |
| **File** | O arquivo. |
| **Fluorochrome** | O nome dele, que vira o nome do canal desmisturado. Obrigatório. |
| **Marker** | O que está marcado com ele, opcional: aparece como o marcador do canal em todo lugar. |
| **Negative** | **Own dim**: os eventos mais fracos do próprio controle. **Unstained**: os eventos do controle não marcado. |

**Tick all** e **Untick all** fazem o que dizem. Abaixo da lista:

| Controle | O que faz |
|---|---|
| **Unstained control** | Células não marcadas ou esferas negativas. É o negativo dos controles postos em **Unstained**, o piso dos pesos no WLS e a fonte do espectro de autofluorescência. Um arquivo com `unstained` ou `negative` no nome é escolhido no começo, de preferência um de células. |
| **Autofluorescence as one more spectrum** | O espectro do próprio não marcado é desmisturado como um canal próprio, chamado `AF`. |
| **Method** | **Ordinary least squares (OLS)** ou **Weighted least squares (WLS)**. |
| **Compute spectra** | Lê os controles e calcula os espectros. |
| **Read CSV…** | Espectros de uma tabela: um fluorocromo por linha e um detector por coluna, ou o contrário (deduzido de onde estão os nomes dos detectores). Cada um é escalado para pico 1. |
| **Save CSV…** | Os espectros, um fluorocromo por linha. |

Os detectores são todos os canais do primeiro controle menos scatter, tempo e semelhantes; desses, só as áreas (`-A`) quando o arquivo tem áreas. O diálogo não deixa escolher os detectores à mão; um CSV de espectros deixa.

Controles na pasta das amostras também são amostras. Desmarque-os na [[samples-tab]] para que fiquem fora dos eventos concatenados.

## Como um espectro é medido

Para cada controle: o detector de pico é aquele em que os eventos mais brilhantes (percentil 99,5) ficam mais acima do negativo. Os positivos são os 5% mais brilhantes nesse detector; os negativos são os 30% mais fracos do mesmo controle, ou os eventos do não marcado. Eventos em 99% da faixa de um detector ou acima ficam de fora, porque a saturação entorta o espectro. O espectro é a mediana dos positivos menos a dos negativos, em cada detector, dividida pelo seu maior valor. Menos de 20 eventos utilizáveis, ou positivos que não são mais brilhantes que os negativos, é um erro que diz qual controle.

O espectro de autofluorescência é a mediana do não marcado em cada detector, dividida pelo seu maior valor.

## Os dois métodos

- **OLS**: `a = y · pinv(M)`, uma matriz para todos os eventos. Rápido.
- **WLS**: cada detector de cada evento recebe o peso 1 / sua variância, tomada como o próprio valor (à maneira de Poisson), nunca abaixo de um piso: o espalhamento robusto (MAD) do não marcado naquele detector, e nunca abaixo de 125. Um detector brilhante, que é ruidoso, passa a contar menos que um fraco. É um sistema pequeno por evento, por isso é lento: cerca de 13 s para 222.000 eventos e 43 espectros, contra 0,02 s do OLS. É o WLS por evento do pacote AutoSpectral, com o mesmo piso de 125.

## As abas

| Aba | O que mostra |
|---|---|
| **Spectra** | Todos os espectros ao longo dos detectores. **Show** destaca um sobre os outros, que ficam apagados. |
| **Similarity** | A similaridade de cosseno de cada par de espectros: 1 é a mesma forma, 0 nada em comum. Pares acima de 0,98 ficam em negrito, porque a Cytek avisa acima disso. Abaixo da tabela, o **índice de complexidade** — o número de condição da matriz de espectros, maior valor singular sobre o menor — e o par mais parecido. Quanto maior o índice, mais o ruído de cada fluorocromo se espalha nos outros. |
| **Check** | Cada controle desmisturado com o conjunto inteiro, nos seus 5% de eventos mais brilhantes: o detector de pico (**Peak**), **Positives**, **Unexplained** (a fração mediana do sinal que os espectros não explicam) e **Spill** com **Into** (a maior mediana que cai noutro canal, como fração do próprio, e em qual canal). Os dois perto de 0 quer dizer que o conjunto explica os seus controles. Clique numa linha para destacar aquele espectro. |

## Aplicar

**Apply** desmistura toda amostra que tem os detectores dos espectros. Os canais desmisturados, com o nome dos fluorocromos, tomam o lugar dos detectores (com suas alturas e larguras); scatter, tempo e o que os espectros não cobrem ficam como estavam. Para o resto do programa eles são fluorescência compensada comum: recebem o eixo biexponencial, `Comp-` nos gates, normalização, estatística, clusters e figuras.

Gates desenhados em detectores que o unmixing substitui não teriam o que cortar, então **Apply** recusa enquanto existirem e diz quais são: apague esses gates antes. **Stop unmixing** volta a compensar as amostras como antes.

Como uma matriz nova, aplicar ou parar muda os valores sob os gates (uma barra na aba Gates pede que sejam revistos) e deixa desatualizados os eventos concatenados. O unmixing é salvo com o projeto.

## Como foi conferido

Medido em 23 de novembro de 2026 num conjunto público do Aurora (42 controles de uma cor em esferas, 64 detectores) que traz cada arquivo bruto e também como o SpectroFlo 3.3 o desmisturou (`docs/comparacao.md`, seção 8):

- Com os espectros do próprio SpectroFlo, o OLS daqui reproduz os valores que o SpectroFlo gravou a menos de 1,3 × 10⁻⁷ do maior valor.
- Os espectros medidos aqui a partir dos mesmos controles têm similaridade de cosseno mediana de 0,99994 com os do SpectroFlo; o pior foi 0,966 (ViaKrome 808).
- OLS, WLS, a similaridade e o número de condição batem com as funções do pacote AutoSpectral a 10⁻⁹ ou melhor.

Não conferido: o WLS e o índice de complexidade da própria Cytek (a fórmula não é publicada, e o conjunto usou OLS), arquivos do Sony ID7000 ou do FACSymphony A5 SE, e amostras de tecido. Não disponível: matriz de espalhamento do spillover, autofluorescência por célula ou vários espectros de autofluorescência, e um unmixing diferente para cada amostra.
