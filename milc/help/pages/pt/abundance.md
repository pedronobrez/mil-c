---
title: Abundância diferencial
---

A abundância diferencial pergunta se um cluster é uma fração maior ou menor da amostra num grupo do que noutro. A pergunta é feita sobre as **contagens** — quantas células de cada amostra caem no cluster — e não sobre as porcentagens, por um modelo binomial negativo à maneira do DA-edgeR do diffcyt (Weber et al., *Commun Biol* 2019).

Abra com **Differential abundance…** na [[clusters-tab|aba Clusters]]. Ele trabalha sobre os clusters como a aba os mostra — juntados, renomeados, com os desenhados à mão — e sobre o agrupamento escolhido lá em **Compare** e **within each**.

## Por que contagens e não porcentagens

Uma porcentagem esquece de onde veio: 0,2% de 109 células e 0,2% de 50.000 são o mesmo número para um teste t, e as porcentagens de um cluster raro estão longe de ser normais. Um modelo de contagens sabe que uma amostra com poucas células dá uma fração menos certa. Em simulações com cinco amostras por grupo ele achou cerca do dobro das mudanças verdadeiras de duas vezes que o teste t de Welch sobre as porcentagens, que a aba Clusters também oferece (ver [[cluster-statistics]]).

## Os controles

| Controle | O que faz |
|---|---|
| **Compare** *coluna* | O grupo do numerador. Uma razão positiva quer dizer mais neste grupo |
| **against** | O grupo de referência |
| **Adjust for** | Opcional: outra coluna das amostras — um lote, um dia, o sexo dos animais. O seu efeito é estimado ao lado do do grupo e retirado, como um termo aditivo do modelo. Só são oferecidas colunas cujo valor muda entre amostras e que não são a de agrupamento nem a de painéis |
| **Test** | **Quasi-likelihood F test** (o padrão) ou **Likelihood-ratio test** |

Toda mudança recalcula na hora. Só as amostras dos dois grupos comparados são usadas, também para a dispersão. Quando a aba Clusters tem uma coluna em **within each**, a análise inteira — modelo, dispersão e tudo — é feita separadamente dentro de cada um dos seus valores, e os p são então ajustados sobre os testes de todos os painéis juntos.

## O método

Para cada cluster, a contagem em cada amostra é modelada como

`contagem ~ NB(média = total × exp(X b), variância = média + dispersão × média²)`

em que *total* é o número de células da amostra em todos os clusters (o offset) e *X* tem um intercepto, o grupo e, se escolhido, o fator ajustado. O coeficiente do grupo é o log da razão, e é ele que se testa.

- **Ajuste** — mínimos quadrados reponderados iterativamente, com o passo reduzido à metade quando piora.
- **Dispersão** — uma dispersão comum a todos os clusters, a que maximiza a soma das suas verossimilhanças de perfil ajustadas de Cox-Reid (o `estimateGLMCommonDisp` do edgeR), procurada numa grade de 2⁻²⁰ a 2³.
- **Quasi-likelihood F test** — o modelo é ajustado na dispersão comum; a quase-dispersão de cada cluster (deviance sobre os graus de liberdade residuais) é puxada em direção às dos outros pelo Bayes empírico do limma, e F = razão de verossimilhança / quase-dispersão encolhida. É o `glmQLFit`/`glmQLFTest` do edgeR anterior à versão 4. Leva em conta que a dispersão é uma estimativa de poucas amostras, e é o mais bem calibrado.
- **Likelihood-ratio test** — o que o `diffcyt::testDA_edgeR` chama: uma dispersão para cada cluster, encolhida em direção à comum com 10 graus de liberdade a priori, e depois um teste qui-quadrado. Acha mais, e é um pouco generoso demais com poucas amostras.
- **Razão** — vem de um segundo ajuste com uma contagem a priori de 0,125 por amostra (proporcional ao seu total), como faz `glmFit(prior.count = 0.125)`. Um cluster ausente de um grupo inteiro fica com uma razão grande mas finita, e uma nota. O teste não usa essa contagem a priori.
- **Ajuste dos p** — Benjamini-Hochberg sobre todos os clusters testados (e todos os painéis).
- Sem fatores de normalização TMM e sem filtro de clusters raros. Um cluster sem nenhuma célula em nenhum dos dois grupos não é testado.

## O resultado

A tabela da esquerda tem uma linha por cluster (e por painel):

| Coluna | Significado |
|---|---|
| *within each* | O painel, quando a aba tem um |
| Cluster | O nome do cluster |
| *Grupo A* (mean %) | A média, sobre as amostras do grupo, da porcentagem da amostra no cluster |
| *Grupo B* (mean %) | O mesmo para a referência |
| log2 fold change | O coeficiente do grupo em log2: 1 é uma fração duas vezes maior, −1 a metade |
| p | O p do teste |
| p adjusted (Benjamini-Hochberg) | A taxa de falsas descobertas |
| (sem nome) | Estrelas pelo p ajustado, como em [[statistics-tab]] |
| n | Amostras em cada grupo, `5 vs 5` |
| Test | O teste usado |
| Statistic | F ou a razão de verossimilhança |
| Dispersion | A do modelo com que o teste foi feito |
| Note | Por que nenhum teste foi feito, ou que um grupo não tinha células |

O **Volcano** à direita desenha cada cluster testado: log2 da razão na horizontal, −log10 do p ajustado na vertical, com uma linha tracejada em FDR 0,05 e linhas-guia em ±1 (duas vezes). Os clusters abaixo de FDR 0,05 são coloridos — "More in *A*" ou "More in *B*" — e nomeados; os outros ficam cinza. Com painéis, um volcano por painel nos mesmos eixos. **Figure…** transforma-o em figura (ver [[figures]]).

Embaixo, um parágrafo diz o que foi feito e o que isso não mostra; ele é parte do resultado. **Export CSV…** (como `differential_abundance.csv`), **Copy** e **Close** ficam no rodapé.

## Limites

- As células contadas são as que foram concatenadas, depois da amostragem, não todas as adquiridas (ver [[concatenation]]).
- O que se testa é a **fração** da amostra que o cluster ocupa: um cluster que cresce torna todos os outros um pouco menores.
- Nas simulações, com cinco amostras por grupo o teste de quase-verossimilhança deu p < 0,05 com a frequência certa quando não havia diferença; com duas ou três por grupo, ou com o teste de razão de verossimilhança, um pouco vezes demais, e valores de p muito pequenos mais ainda. Duas amostras por grupo não bastam para nenhum dos dois testes.
- Posto ao lado do edgeR 4.10 chamado do mesmo jeito, dá os mesmos números; contra os padrões atuais do edgeR e do diffcyt, de 10% a 20% das descobertas diferem (ver [[measured-facts]]).

Para saber se as células de um cluster expressam mais ou menos de um marcador, ver [[differential-expression]].
