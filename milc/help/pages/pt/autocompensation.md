---
title: Autocompensação por controles de uma cor
---

Uma matriz de spillover pode ser calculada a partir de controles de uma cor: um arquivo FCS por fluorocromo, cada um marcado só com aquele fluorocromo, e, se houver, um não marcado. Num controle marcado só com o fluorocromo *i*, cada detector é uma reta contra o detector do próprio controle, o primário; a inclinação dessa reta é o spillover de *i* naquele detector, que é a linha *i* da matriz (veja [[compensation]]).

## Onde fica

Na [[gates-tab]], **Compensation…** e depois **From controls…**. O diálogo **Compensation from single-stain controls** calcula uma matriz para os detectores da matriz mostrada no diálogo Compensation. Calcular não aplica nada: o resultado vai para a tabela de compensação, de onde é dado às amostras escolhidas, como uma matriz digitada à mão.

O mesmo existe na linha de comando: `milc autocomp CONTROLES --out spillover.csv` (veja [[command-line]]).

## Os controles

| Controle | O que faz |
|---|---|
| **Add files…** | Escolhe os arquivos FCS dos controles. São lidos fora da linha de execução da janela. |
| **Remove** | Tira da lista os controles selecionados. |
| Método | **AutoSpill (robust, refined)** ou **Classic (positive vs negative)**, descritos abaixo. |
| **Gate on scatter** | Fica com a população principal em forward contra side scatter, como a elipse robusta que contém 95% dela ([[auto-gates]]): debris e grumos ficam de fora. |
| **Unstained as negative** | Só no clássico: os eventos do controle não marcado são o negativo de todos os controles, no lugar dos eventos negativos de cada controle. |
| **Compute** | Calcula a matriz a partir dos controles da lista. |

A tabela lista cada **Control**, seus **Events**, seu **Detector** e o detector em que ele é mais brilhante (**Brightest in**). Cada controle recebe o detector em que é mais brilhante, e nenhum detector é dado duas vezes; troque na lista **Detector** onde estiver errado, ou escolha **Unstained**. Um arquivo cujo nome diz que é não marcado (`unstained`, `no stain`, `blank`, `negative control`, ou `US` ou `neg` como palavra) começa como o não marcado. Usam-se até 100.000 eventos de cada controle, sorteados.

O brilho é quanto os eventos mais brilhantes de um controle (o percentil 99,5) ficam acima da mediana do controle não marcado em cada detector, ou acima do percentil 20 do próprio controle quando não há não marcado. Os detectores são então casados com os controles de modo que cada controle receba aquele em que é relativamente mais brilhante.

## O resultado

A linha de status diz quantos detectores foram calculados e, no AutoSpill, quantas iterações rodaram, se convergiram, e o maior spillover que sobra nos controles depois de compensados com a matriz achada. As notas abaixo listam o que saiu diferente (um ajuste robusto que não convergiu, um controle sem gate de scatter e por quê, detectores sem controle).

Duas matrizes aparecem lado a lado, em porcentagem: **From the controls** e a matriz que as amostras têm (com o título do lugar de onde ela vem: o arquivo da amostra ou o workspace do FlowJo). Uma célula colorida difere em mais de 1 ponto percentual entre as duas; a dica dela dá o valor das amostras. As linhas dos detectores sem controle são copiadas da matriz das amostras.

| Botão | O que faz |
|---|---|
| **Put in the compensation table** | Fecha este diálogo e digita a matriz achada na tabela de compensação, com quatro algarismos significativos. Depois selecione as amostras e aplique lá. |
| **Save as CSV…** | A matriz em proporções (não porcentagens), com os detectores na primeira linha e na primeira coluna, como a keyword `$SPILLOVER` guarda. |
| **Cancel** | Fecha sem mexer na tabela. |

## Os dois métodos

### Clássico

O que os programas dos citômetros sempre fizeram. Em cada controle, o vale entre os eventos negativos e positivos do detector primário é achado na escala biexponencial (duas médias, depois o ponto mais baixo de um histograma suavizado entre elas; um vale mais alto que 80% do pico menor não é vale). Os positivos são os que ficam acima dele, os negativos os que ficam abaixo — as mesmas células, o mesmo tubo, a mesma autofluorescência. O 1% de cima e o de baixo ficam de fora. A linha é a mediana dos positivos menos a dos negativos, em cada detector, dividida pela mesma diferença no detector primário.

Os eventos do não marcado são usados como negativo quando **Unstained as negative** está marcado, ou quando um controle tem menos de 50 eventos negativos próprios (esferas muitas vezes não têm nenhum). Um controle com uma população só, havendo um não marcado, é dividido no percentil 99,5 do não marcado. São necessários pelo menos 50 eventos positivos e 50 negativos.

### AutoSpill

O método de Roca et al. (2021, *Nature Communications* 12: 2890), escrito a partir do pacote em R dos autores (`autospill` 0.2.0).

1. Uma regressão linear robusta (o estimador M de Huber, como o `MASS::rlm`) de cada detector contra o primário, em todos os eventos do controle, deixando de fora o 1% de cima e o de baixo, dá uma primeira matriz. Onde o ajuste robusto não converge, usam-se mínimos quadrados comuns.
2. Os controles são compensados com ela, a mesma regressão nos eventos compensados mede o spillover que sobrou, e a matriz é corrigida por ele — de novo e de novo. Primeiro na escala linear, até nenhuma inclinação restante passar de 0,01; depois na biexponencial, até nenhuma passar de 0,0001, ou até nada mais melhorar. O passo é a correção inteira no começo e um décimo dela quando o progresso para; no máximo 100 iterações.

Os limiares, os passos e as regras de parada do pacote foram mantidos como são.

## Como foi conferido

Medido em 23 de novembro de 2026 contra o AutoSpill 0.2.0 e o flowStats 4.24.0 (`docs/comparacao.md`, seção 5):

- Nos mesmos eventos, o AutoSpill daqui dá a matriz do pacote a menos de 2 × 10⁻⁵ em controles simulados e de 3 × 10⁻⁴ nos controles reais que vêm com o flowCore (um FACSCalibur, quatro detectores), com o mesmo número de iterações.
- Em controles simulados com matriz conhecida, acha a matriz verdadeira a menos de 0,04 ponto percentual; o método clássico, a menos de 0,3.
- Com 15% de debris nos controles, o gate de scatter levou o erro do AutoSpill de 0,26 para 0,02 ponto.

Não conferido: controles de um citômetro moderno de vários lasers, esferas, citômetros espectrais, os controles do artigo do AutoSpill (não puderam ser baixados), e nenhuma matriz calculada aqui foi posta ao lado da do FlowJo ou do FACSDiva para os mesmos controles. Antes de confiar numa matriz para um experimento, compare-a com a que o programa de aquisição calculou e confira o detector que cada controle recebeu.
