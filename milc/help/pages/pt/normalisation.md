---
title: Normalização de lotes
---

Amostras adquiridas em dias diferentes diferem por razões que não são biologia: um laser, um lote
de anticorpo, uma voltagem, uma limpeza. Se cada lote tiver amostras de referência que deveriam
ser iguais — alíquotas de um mesmo controle, marcadas e adquiridas com cada lote — então o que
separa as referências entre os lotes é técnico, e pode ser aprendido e retirado de todas as
amostras do lote. O MIL-C faz isso à maneira do CytoNorm (Van Gassen et al., Cytometry A 2020),
antes de as amostras serem concatenadas.

## Os ajustes

Na seção **Batch normalisation** da [[samples-tab]]. As colunas oferecidas são as da tabela da
aba Samples, então primeiro acrescente uma coluna que diga em que lote cada amostra foi
adquirida.

| Controle | Significado |
|---|---|
| **Batches** | a coluna que diz o lote de cada amostra — um dia, uma corrida. **None** desliga a normalização |
| **Reference** | uma coluna e um dos seus valores: as amostras que deveriam ser iguais em todo lote, como Group = Control |
| **Within** | `each marker whole` (1, cada marcador inteiro), ou um número de clusters do FlowSOM (2 a 30) dentro dos quais a correção é ajustada |
| **Bring to** | `The mean of the batches` (a média dos lotes), ou um lote, `Day 1, left as acquired`, para o qual os outros são trazidos |
| **Also in Gates and Statistics** | normaliza as amostras ao serem lidas, não só os eventos concatenados |

Os outros controles ficam ativos quando **Batches** indica uma coluna.

Escolha a referência com cuidado. O que separa as amostras de referência entre os lotes é retirado
de todas as amostras do lote — inclusive qualquer diferença real entre elas. Alíquotas de uma
mesma amostra controle são a escolha segura; um grupo controle de animais diferentes é uma
suposição de que o grupo não muda entre os lotes.

## Antes e depois de concatenar

Abaixo dos ajustes, uma linha diz o que será feito:

> At the next Concatenate: 1: 2 reference sample(s); 2: 2 reference sample(s). Whatever
> separates the Control samples between batches will be removed from every sample — assuming
> they should not differ.

isto é, quantas amostras de referência cada lote tem, e o aviso de que tudo o que separa as
referências será retirado. Ou diz o que atrapalha:

- `Needs samples with Group = Control in at least two batches; there are such samples in 1.`
  (são precisas referências em pelo menos dois lotes)
- `Batch 3, which the others are to be brought onto, has no sample with Group = Control.` (o lote
  alvo não tem referência)

Uma amostra de referência só conta se estiver marcada e tiver a população escolhida.

Depois de **Concatenate**, a linha relata o resultado, por exemplo:

> Normalised. References apart by 0.0220 before and 0.0135 after, on average over the markers
> (share of the axis, on reference events held out). Most moved: CD4, 0.0410 → 0.0120.

Ou seja, a distância média entre as referências antes e depois, como fração do eixo, medida em
eventos separados do ajuste, e o marcador que mais se moveu. A linha diz ainda quantos pares de
lote e cluster tinham poucos eventos de referência e ficaram como estavam, e quais amostras não
foram normalizadas porque o seu lote não tinha referência.

## Olhando: Before and after…

**Before and after…** aparece depois de uma normalização. Ele abre **Batch normalisation: before
and after**: para o **Marker** escolhido, os eventos de referência de cada lote desenhados uns
sobre os outros como histogramas suavizados, antes (em cima) e depois (embaixo), cada curva na
escala do seu próprio pico. Onde as curvas se sobrepõem, os lotes concordam.

Os eventos mostrados são separados do ajuste: metade dos eventos de referência de cada lote é
posta de lado, as curvas são ajustadas na outra metade, e o "depois" é medido na metade que elas
nunca viram. Medida nos eventos em que foi ajustada, a concordância seria perfeita por construção
e não diria nada. **Figure: before…** e **Figure: after…** transformam qualquer dos dois painéis
em figura ([[figures]]).

## O que é feito

Para cada lote e cada marcador:

1. Os eventos de referência do lote são reunidos — todos os eventos da população escolhida de
   cada amostra de referência daquele lote, lidos dos arquivos, não só os sorteados.
2. Os seus 99 quantis, de 1% a 99%, são tomados na escala de exibição do marcador (a escolhida em
   [[scales-and-transforms]], senão a biexponencial do FlowJo).
3. O alvo é a média dos quantis dos lotes, ou os do lote escolhido.
4. Uma curva monótona (um spline de Fritsch-Carlson) pelos pares (quantil do lote, quantil do
   alvo) é aplicada a todo evento de toda amostra do lote. Ela nunca inverte a ordem de dois
   eventos.

Só os canais compensados — os marcadores de fluorescência — são normalizados; dispersão e tempo
ficam como estão.

Com **Within** acima de 1, um mapa do FlowSOM (10 × 10, treinado com o mesmo número de eventos de
referência de cada lote) é cortado em tantos clusters, as curvas são ajustadas dentro de cada
cluster, e cada evento é corrigido pela curva do cluster em que cai. Isso serve para um desvio que
difere entre tipos celulares. Exige muitos eventos de referência: um lote com menos de 50 num
cluster não é corrigido ali e aparece listado no resultado.

Além dos quantis extremos nada se sabe sobre o efeito de lote, então a correção continua como o
deslocamento que era no último quantil, e não pela inclinação final do spline.

## Também em Gates e Statistics

Por padrão só os eventos concatenados são normalizados, e as abas Gates e Statistics mostram os
valores como adquiridos. Marcar **Also in Gates and Statistics** normaliza cada amostra ao ser
lida: gráficos, gates, contagens e a tabela de estatísticas passam a usar valores normalizados.
Gates desenhados antes foram desenhados sobre os valores como adquiridos, então a aba Gates pede
que sejam revistos ([[gates-tab]]), e os eventos concatenados ficam marcados como desatualizados
até você concatenar de novo. As curvas são sempre aprendidas a partir dos arquivos como
adquiridos, e reajustadas a cada **Concatenate**, não enquanto você edita gates: um gráfico não
se mexe debaixo de um gate enquanto ele é desenhado.

## Comparado com o CytoNorm

Posto ao lado do pacote em R (CytoNorm 2.0.13) em 23 de novembro de 2026, nos mesmos eventos:

- Sem clusters, as duas correções diferem, dentro do 1%–99% das referências, por 1e-5 a 2e-5 do
  eixo em média, onde a correção média é de 0,02 a 0,04: é a mesma curva.
- Com dez clusters os mapas dos dois programas não são o mesmo mapa, e os resultados ficam tão
  distantes quanto duas sementes de qualquer um deles; os dois deixam as referências igualmente
  próximas (0,013 a 0,014 do eixo, partindo de 0,022).
- Fora da faixa das referências eles se separam de propósito. A inclinação final do CytoNorm
  moveu valores isolados em até 1,1 de um eixo que vai de 0 a 1; a do MIL-C, em 0,17 no máximo.

Os detalhes estão em [[measured-facts]].

## Pela linha de comando

```
milc concat DADOS --wsp analysis.wsp --normalize-by Day --reference Group=Control \
    --normalize-clusters 1 --normalize-to 1 --out concat
```

`--normalize-to` é opcional (padrão: a média dos lotes). O comando imprime a distância entre as
referências antes e depois. Veja [[command-line]].
