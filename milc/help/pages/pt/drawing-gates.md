---
title: Desenhar e editar gates
---

Todo gate é desenhado na população mostrada na [[gates-tab]], sobre os dois parâmetros do gráfico (ou sobre o único parâmetro de um histograma). Ele vira filho dessa população na árvore. Escolha primeiro a população, depois os eixos, depois a ferramenta.

## As ferramentas

Os nove botões sob **New gate on the population shown**:

| Botão | O que faz |
|---|---|
| **Rectangle** | Arraste sobre os eventos. Num histograma, desenha uma faixa de X. |
| **Polygon** | Clique em cada vértice; clique duas vezes para fechar (no mínimo três vértices). |
| **Lasso** | Segure o botão do mouse e contorne os eventos à mão livre; solte para fechar. Vira um polígono de algumas dezenas de vértices, que se move como qualquer outro. |
| **Quadrant** | Clique onde as duas linhas devem se cruzar: quatro populações de uma vez. |
| **Ellipse** | Arraste a caixa em que a elipse cabe; depois gire e estique pelas alças. |
| **Auto…** | Um menu de gates propostos a partir dos eventos. Veja [[auto-gates]]. |
| **Combine…** | Uma população feita de outras da população mostrada: em todas, em alguma, ou em nenhuma. |
| **Rename** | Um nome novo para a população selecionada na árvore. |
| **Delete** | Apaga a população selecionada e tudo o que está dentro dela. |

As ferramentas de desenho ficam apertadas até o gate ser desenhado; apertar outra ferramenta, apertar a mesma de novo ou **Esc** cancela. Toda ferramenta, menos **Rectangle**, precisa de dois parâmetros: num histograma, a barra de status diz `That gate needs two parameters: choose a Y axis`.

Quando a forma se fecha, uma caixa pede o nome. O nome sugerido é `Gate`, ou `Gate 2`, `Gate 3`… quando já existe: duas populações do mesmo pai não podem ter o mesmo nome.

## O que cada gate guarda

- Um **retângulo** fica nas unidades em que o arquivo guarda os valores, então não depende da escala dos eixos. Um evento na borda de baixo está dentro; um na borda de cima, fora. Uma faixa num histograma é um retângulo em um parâmetro.
- Um **polígono** (e um laço) é testado na escala em que foi desenhado: um lado reto num eixo biexponencial é a linha com que os eventos são comparados, como faz o FlowJo. Se depois a escala daquele eixo mudar, o polígono acompanha: os vértices ficam nos mesmos valores e os lados ficam retos na escala nova. Veja [[scales-and-transforms]].
- Uma **elipse** guarda centro, semieixos e ângulo nas coordenadas em que foi desenhada. Uma elipse não pode ir para outra escala sem virar outra forma; então, depois que o eixo dela recebe outra escala, ela aparece como a curva que virou, com nome e fração, mas não pode mais ser arrastada.
- Um **quadrante** corta o gráfico num ponto, nas unidades do arquivo. Um evento exatamente sobre uma linha vai para o lado alto.

## O laço

O que a mão desenha tem centenas de pontos. Ele é reduzido a um polígono pelo método de Ramer–Douglas–Peucker: só ficam os vértices sem os quais o traçado se afastaria do polígono mais que 0,4% da largura ou da altura do gráfico. Se ainda sobrarem mais de 60 vértices, a tolerância é afrouxada em 50% de cada vez, até haver 60 ou menos. O mouse parado num lugar não acrescenta nada. Onde a mão voltou sobre o próprio traço, o polígono acompanha.

## Quadrantes

Um quadrante cria quatro populações, com o nome feito da primeira palavra do rótulo de cada eixo e do lado de cada linha em que estão: `CD4+ CD8-`, `CD4+ CD8+`, `CD4- CD8+`, `CD4- CD8-`. Internamente são numerados como o FlowJo numera: Q1 em cima à esquerda, Q2 em cima à direita, Q3 embaixo à direita, Q4 embaixo à esquerda. O nome e a fração de cada canto aparecem no seu canto do gráfico.

- Arraste qualquer das duas linhas para mover o cruzamento; as quatro populações se movem juntas.
- Cada braço tem uma alça quadrada pequena. Arraste uma para inclinar só aquele braço. A partir daí os quatro braços são semirretas independentes (o quadrante do FlowJo com braços deslocados), e uma alça redonda no centro move os quatro.
- Um braço não pode ser girado até o braço vizinho nem além dele; se for solto ali, volta, e a barra de status avisa. Todo evento está sempre em exatamente um canto.

O **quadrante curly** do FlowJo, cujos braços de cima e da direita se curvam para longe do centro, é lido de um workspace mas não é desenhado aqui; pode ser movido pelo centro. Veja [[workspaces-and-gating-ml]].

## Mover, renomear, apagar

- Arraste o contorno de um gate para movê-lo, ou as alças para mudar a forma. As contagens na árvore e as porcentagens no gráfico se atualizam quando o mouse é solto.
- Com **Apply to every sample** marcado, mover, renomear ou apagar vale para todas as amostras. Desmarcado, só para a mostrada.
- **Rename** leva o nome novo a todas as amostras (quando marcado), e uma combinação feita da população é renomeada junto.
- **Delete** remove a população, todas as populações abaixo dela e toda combinação ao lado que era feita dela (que de outro modo ficaria vazia). Apagar a população criada pelo controle de qualidade desliga o controle de qualidade (veja [[quality-control]]).
- **Back to the workspace's gates** (sob as ferramentas) esquece tudo o que foi desenhado ou movido e volta aos gates do workspace.

## Combinações

**Combine…** abre **Combine populations**. Marque as populações da mostrada que serão combinadas e escolha em **Events**:

| Escolha | Eventos que ficam |
|---|---|
| **In all of them (AND)** | dentro de todas as populações marcadas |
| **In any of them (OR)** | dentro de pelo menos uma |
| **In none of them (NOT)** | dentro de nenhuma |

A caixa **Name** sugere um nome feito da escolha, como `CD4+ and CD69+` ou `not CD4+ or CD8+`. Uma combinação não tem forma e não aparece no gráfico. Se uma população que ela cita deixar de existir, a combinação fica vazia em vez de adivinhar.

## Desenhar em várias amostras

Selecionar várias amostras em **Samples** sobrepõe todas no gráfico, mas os gates mostrados e desenhados são os da primeira amostra (a legenda diz isso). Com **Apply to every sample** marcado, isso não faz diferença: todas as amostras recebem o gate.

Para pôr um gate num lugar diferente em cada amostra, desmarque **Apply to every sample** e desenhe ou mova o gate em cada amostra, ou use **Fit on each sample separately** no menu de [[auto-gates]].
