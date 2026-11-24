---
title: Escalas e transformações
---

Uma escala é o modo como um valor medido vira uma posição num eixo. A fluorescência compensada tem valores negativos e se estende por cinco décadas: um eixo linear esmaga toda população fraca contra o zero, e um logarítmico não consegue mostrar os negativos. As escalas abaixo são as respostas de costume. No MIL-C a escala pertence ao parâmetro, não ao gráfico: todo gráfico daquele parâmetro, em todas as abas, usa a mesma, e assim também os algoritmos de [[embeddings]] e [[clustering-methods]], que recebem cada marcador na sua escala.

## Em que escala um parâmetro começa

- A fluorescência (todo canal que a matriz de compensação cobre, ou todo canal desmisturado) começa na biexponencial do FlowJo com os ajustes de costume.
- Scatter, tempo e qualquer outro parâmetro começam lineares, como guardados.

## Como mudar

Aperte **Scale…** ao lado de **X** ou de **Y** na [[gates-tab]]. O diálogo **Scale of** *parâmetro* tem:

| Controle | O que faz |
|---|---|
| **Scale** | **Biexponential**, **Logicle**, **Arcsinh**, **Logarithmic** ou **Linear**. Os ajustes abaixo mudam com ela. |
| **Every fluorescence parameter, not just this one** | Dá a mesma escala a todos os parâmetros de fluorescência de uma vez. |
| **The usual one** | Volta à escala inicial descrita acima, e fecha. |
| **OK** / **Cancel** | Aplica a escala, ou deixa como estava. |

As escalas são salvas com o projeto (veja [[projects-and-saving]]).

## As escalas e seus ajustes

### Biexponencial

A biexponencial do FlowJo, calculada como o cytolib traduziu a rotina do FlowJo: uma tabela de 4.096 pontos lida por interpolação, de 0 na base do eixo a 1 no topo, com o topo em 262.144.

| Ajuste | Faixa | Padrão | Significado |
|---|---|---|---|
| **Linear width** | 1 a 5.000 | 100 | A largura da parte quase linear em torno do zero. Maior junta uma população negativa espalhada. |
| **Negative decades** | 0 a 3 | 0 | Mais espaço abaixo do zero, em décadas. |
| **Positive decades** | 2 a 7 | 4,4185 | As décadas do lado positivo, o padrão do FlowJo. |

A forma dela (monotônica, zero onde deve, a faixa toda no topo) é testada. Ela **não** foi comparada com valores tirados do próprio FlowJo.

### Logicle

A Logicle publicada (Parks, Roederer e Moore, 2006), nos parâmetros de Moore e Parks (2012). Conferida com a tradução em C do código dos próprios autores feita pelo FlowUtils.

| Ajuste | Padrão | Significado |
|---|---|---|
| **Linear width (W), decades** | 0,5 | As décadas dadas à parte quase linear em torno do zero. Entre 0 e M / 2. |
| **Decades (M)** | 4,5 | As décadas que a escala inteira abrange. Entre 1 e 8. |
| **Negative decades (A)** | 0 | Mais espaço abaixo do zero. Entre −W e M − 2W. |
| **Top of scale (T)** | 262.144 | O valor no topo do eixo; 262.144 para um citômetro de 18 bits. |

O diálogo mantém W e A dentro das faixas permitidas conforme M e W mudam.

### Arcsinh

`asinh(x / cofator)`. **Cofactor** de 1 a 10.000, padrão 150, o valor de costume para fluorescência convencional; 5 é o de costume para citometria de massa. O eixo vai do valor de −10 × cofator, para dar espaço a uma população negativa compensada, até o valor de 262.144.

### Logarítmica

Um logaritmo simples de **Lowest value** (0,01 a 10.000, padrão 10) até 262.144. O zero e todo valor negativo ficam na borda de baixo — que é a razão de existirem as escalas biexponenciais.

### Linear

Os valores como guardados, de zero até a faixa do parâmetro no arquivo.

Em toda escala que não a linear, o eixo é marcado em 0, 10², 10³… onde há espaço para o rótulo.

## O que muda quando uma escala muda

- **Polígonos**, laços e polígonos propostos por **Auto…** acompanham o eixo: os vértices ficam nos mesmos valores e os lados ficam retos na escala nova, que é onde passam a ser testados. As contagens podem mudar um pouco por isso.
- **Braços de quadrante** inclinados acompanham do mesmo jeito; um quadrante curly se curva nos canais da escala nova.
- **Retângulos**, faixas e quadrantes simples ficam nas unidades do arquivo e não se movem.
- **Elipses** ficam na escala em que foram desenhadas. Aparecem como a curva que viraram no eixo novo e não podem mais ser arrastadas; volte à escala antiga para mover uma.
- Os gates ajustados por **Auto…** são ajustados de novo nos valores novos (veja [[auto-gates]]).
- Todas as contagens são refeitas, e as outras abas redesenham os gráficos daquele parâmetro.

A nota no pé do diálogo diz o mesmo: os polígonos desenhados neste parâmetro são testados na escala nova, e os algoritmos de Explore recebem cada marcador na sua escala.

## Qual usar

Para a maioria dos painéis convencionais e espectrais, a biexponencial padrão é o começo certo. Aumente **Linear width** (ou W) quando uma população negativa estiver espalhada numa faixa larga, e diminua quando duas populações fracas se misturarem. Arcsinh com cofator 150 é o que muitos fluxos em R usam, e dá figuras parecidas. Use a logarítmica só para parâmetros que nunca são negativos. Scatter fica melhor linear. A compensação em si está em [[compensation]].
