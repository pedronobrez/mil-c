---
title: Gates propostos a partir dos eventos
---

**Auto…** na [[gates-tab]] propõe um gate a partir dos eventos da população mostrada: o vale entre dois picos de um histograma, o contorno ou a elipse de uma população em dois parâmetros, ou a faixa de singletos ao longo de uma diagonal. A proposta aparece tracejada no gráfico primeiro. Aceita, é uma faixa, um polígono ou uma elipse comum: pode ser movida, renomeada e apagada como uma desenhada à mão (veja [[drawing-gates]]). Quando os eventos não mostram o que o método procura, nada é proposto, e uma mensagem diz por quê.

Tudo é medido nas coordenadas do próprio gráfico, depois da escala de cada eixo (veja [[scales-and-transforms]]): uma população é um pico onde ela é olhada, não nas unidades em que o arquivo guarda os valores.

## O menu

| Item | O que é proposto |
|---|---|
| **Valley on X: the positive side** | Uma faixa de X do vale para cima. |
| **Valley on X: the negative side** | Uma faixa de X do vale para baixo. |
| **Population at a click: contour** | Clique numa população: o contorno que contém 90% dos seus eventos, como polígono. |
| **Population at a click: ellipse** | Clique numa população: a elipse de uma gaussiana ajustada ao seu núcleo. |
| **Largest population: contour** | O mesmo contorno, para a população com mais eventos. |
| **Largest population: ellipse** | A mesma elipse, para a população com mais eventos. |
| **Singlets: band along the diagonal** | Uma faixa de três desvios-padrão de cada lado de uma reta de Y contra X. |
| **Singlets: a tighter band** | O mesmo, com o espalhamento medido só do lado mais limpo da reta. |
| **Fit on each sample separately** | Uma opção, marcada ou não, para o próximo gate aceito: veja abaixo. |

Depois da proposta, uma caixa mostra o que foi achado (por exemplo `valley at 18% of the lower peak; 32.4% of the events above it`) e pede um nome. O nome sugerido é a primeira palavra do rótulo de X com `+` ou `-` para um vale, `Population` para contorno ou elipse, e `Singlets` para a faixa. **OK** aceita; **Cancel** descarta a proposta.

O vale usa só X e funciona num histograma. Os outros precisam de dois parâmetros.

## Como cada um é achado

### Vale entre dois picos

1. O histograma de X é feito onde os eventos estão: do ponto de 0,1% ao de 99,9% deles, alargado em 5%, em 256 classes. Os eventos empilhados nas pontas do eixo ficam de fora.
2. É suavizado por uma gaussiana cuja largura depende do número de eventos: 3 classes para 5.000 eventos ou mais, até 10 classes para poucos eventos.
3. Um pico só conta se se destaca dos arredores por cinco vezes o que o ruído de contagem sozinho faz numa classe suavizada (e por pelo menos 0,1% do mais alto). Uma população pequena é um pico baixo e largo, e conta assim mesmo.
4. Os dois picos que mais se destacam são as duas populações; o limiar é o ponto mais baixo entre eles.

Nada é proposto quando há menos de 100 eventos, quando há um pico só, quando o vale ainda está acima de 70% do pico mais baixo, ou quando um dos lados ficaria com menos de 1% dos eventos. O limiar, então, é questão de julgamento, ou de um controle.

### Populações em dois parâmetros

Primeiro o gráfico é dividido em populações. Os eventos são contados numa grade de 128 × 128 (as células das bordas, onde os eventos se empilham, são zeradas) e suavizados. Cada célula é ligada ao pico a que se chega subindo a partir dela. Um pico cuja passagem para um vizinho é mais alta que 85% da própria altura, ou cujo terreno tem menos de 1% dos eventos, é uma ondulação do vizinho e se junta a ele. Cada pico que sobra, com seu terreno, é uma população.

- **Contorno**: as células da população escolhida que contêm 90% dos seus eventos, das mais densas para as menos, e a linha em volta delas. Como o gráfico foi dividido antes, a linha para no vale até o vizinho em vez de contornar os dois. É reduzida a no máximo 32 vértices, para poder ser movida à mão. A mensagem diz quantos eventos ficaram dentro e quantas outras populações há no gráfico.
- **Elipse**: uma gaussiana bidimensional ajustada ao núcleo da população. O ajuste começa da metade mais densa da população e repete: fica com os eventos dentro da distância que contém 90% de uma gaussiana, toma média e covariância (aumentada pelo que o corte das caudas tira), até nada mudar. Eventos soltos e outras populações ficam fora do corte e não puxam a elipse. A elipse desenhada contém 95% dessa gaussiana. É para isso que se usam, no R, o `flowCore::norm2Filter` ou um `flowClust` de um cluster.

Um clique que cai entre populações escolhe o pico mais próximo.

### Faixa de singletos

Para singletos na área de um pulso contra a altura: um dubleto tem mais área para a sua altura e sai da reta. A reta é ajustada pela regressão biweight de Tukey, que não dá voz aos eventos longe dela. O espalhamento em volta dela é a mediana das distâncias absolutas dos eventos que o ajuste considera, convertida em desvio-padrão. A faixa vai três desses para cada lado, do menor ao maior X dos eventos que contém. É o método do `flowStats::singletGate` (um modelo linear robusto e uma faixa de predição).

**A tighter band** serve para uma reta com um leque de dubletos de um lado, que alarga o espalhamento: a reta é levada para a crista, onde os eventos são mais densos de través, e o espalhamento é medido só do lado mais limpo.

Nada é proposto quando os eventos não seguem uma reta: quando o espalhamento em volta da melhor reta passa da metade do espalhamento do próprio Y.

## Um gate para todas, ou um para cada amostra

Com **Apply to every sample** marcado e **Fit on each sample separately** desmarcado, o gate é ajustado na amostra mostrada e o mesmo gate vai para todas. Fica registrado que foi ajustado ali: quando a população em que ele está muda (um gate acima foi movido, uma compensação nova, uma escala nova), ele é ajustado de novo naquela amostra e dado a todas outra vez.

Com os dois marcados, cada amostra recebe o gate ajustado nos seus próprios eventos, e ele fica onde a população daquela amostra está. Cada amostra passa a ter o seu: movê-lo numa amostra move só ali. Onde nada pode ser achado numa amostra, ela fica com o gate da amostra mostrada, e a barra de status e a dica da população dizem isso (`not fitted to this sample: …`).

Um gate ajustado que foi movido à mão não é ajustado de novo sozinho quando a população muda, para não perder o trabalho manual sem perguntar. A barra acima do gráfico diz quais gates estão nesse estado; **Refit** ajusta de novo, **Keep** deixa onde estão.

## Limites

São propostas, não substituem olhar. Um marcador sem vale não recebe gate de vale; uma população no ombro de outra pode ser juntada a ela; um contorno é um contorno de densidade, não de biologia. Olhe cada proposta antes de aceitar, e os gates de cada amostra depois de **Fit on each sample separately**. A mesma elipse robusta é o gate de scatter da [[autocompensation]].
