---
title: Comparar clusters entre grupos
---

Um cluster é uma fração maior das amostras num grupo do que noutro, mais do que o acaso explicaria? A [[clusters-tab|aba Clusters]] responde de três maneiras:

- o **Test** na segunda linha da aba, sobre a porcentagem de cada amostra em cada cluster — descrito aqui;
- **Differential abundance…**, sobre as contagens, por um modelo binomial negativo: ver [[abundance]];
- **Differential expression…**, sobre o que as células de cada cluster expressam: ver [[differential-expression]].

Os mesmos testes, sobre qualquer coluna da tabela de estatísticas, estão em **Compare groups…** na [[statistics-tab|aba Statistics]].

## A unidade é a amostra

Todo teste aqui toma um valor por amostra — a fração dos eventos da amostra no cluster — e compara esses valores entre os grupos. Cinco animais são cinco observações, quantas células cada um tenha dado. Testar células contra células acharia "significativa" qualquer diferença e não diria nada sobre os animais. O n de um grupo é o seu número de amostras.

A fração é dos eventos **concatenados** da amostra: o que a [[concatenation|concatenação]] tirou dela, depois da amostragem.

## Os testes

Escolha em **Test**:

| Opção | Dois grupos | Três ou mais |
|---|---|---|
| **Welch t test / ANOVA** | Teste t de Welch, que não supõe variâncias iguais | ANOVA de um fator |
| **Mann-Whitney / Kruskal-Wallis** | Teste U de Mann-Whitney (bilateral) | Kruskal-Wallis |
| **Welch t test / Welch ANOVA** | Teste t de Welch | ANOVA de Welch (Welch, 1951): cada grupo pesa n / variância, para grupos com variâncias diferentes |

O primeiro é paramétrico e o mais poderoso quando os valores são aproximadamente normais dentro de cada grupo; os testes de postos não supõem nada sobre a forma, mas não conseguem dar um p pequeno com poucas amostras. Com três amostras contra três, o Mann-Whitney nunca fica abaixo de p = 0,1, sejam quais forem os dados; a Note avisa ("smallest p possible with these n").

Com mais de dois grupos, **Results…** (e **Compare groups…** em Statistics, quando a sua caixa está marcada) testa também cada par de grupos, ao lado do teste de todos juntos. O gráfico mostra só o teste de todos os grupos.

### Desenhos pareados

Quando os mesmos animais ou doadores foram medidos em todos os grupos — antes e depois, ou várias condições do mesmo doador — escolha em **Paired by** a coluna que diz de qual deles é cada amostra. Os grupos passam a ser condições e o teste é da mudança de cada sujeito:

| Opção | Dois grupos | Três ou mais |
|---|---|---|
| Welch t (qualquer das duas) | Teste t pareado | ANOVA de medidas repetidas, corrigida por Greenhouse-Geisser |
| Mann-Whitney / Kruskal-Wallis | Teste dos postos sinalizados de Wilcoxon | Teste de Friedman |

Só entram sujeitos com valor em todos os grupos comparados; os outros ficam de fora e são contados na Note. Se um sujeito aparece duas vezes num grupo, ou nenhum sujeito está em todos os grupos, não há teste, e um aviso ao lado dos controles diz por quê: "No test: the samples cannot be paired by Donor — …".

A ANOVA de medidas repetidas retira o que difere entre sujeitos antes de comparar os grupos. O seu teste F supõe que as diferenças entre quaisquer dois grupos se espalham igualmente (esfericidade) e dá um p pequeno demais quando não é assim; o ε de Greenhouse-Geisser (de 1, quando isso vale, até 1/(k − 1)) multiplica os dois graus de liberdade. A Note dá o ε e o p sem a correção.

## Ajustar por muitos testes

Cada cluster é um teste, e com vinte clusters um sairá "significativo" a 0,05 só por acaso. Os valores de p são por isso ajustados sobre todos os testes feitos juntos — cada cluster, cada painel, cada par:

| Opção | O que controla |
|---|---|
| **Benjamini-Hochberg** (padrão) | A taxa de falsas descobertas: entre os clusters chamados significativos, a fração esperada que não é. O p ajustado costuma ser chamado de FDR ou valor q |
| **Holm** | O erro por família: a chance de qualquer falso positivo. Mais estrito que o Benjamini-Hochberg, nunca mais estrito que o Bonferroni |
| **Bonferroni** | O mesmo, de modo mais simples: p × o número de testes, no máximo 1 |
| **None** | Sem ajuste |

Testes que não puderam ser feitos não contam como testes.

## No gráfico

Com um teste escolhido, as barras de cada cluster no painel **Frequency** levam o p ajustado do teste de todos os grupos — `p = 0.003`, `p < 0.0001` — onde ele fica abaixo de 0,05, e `ns` onde não fica. Todos os números estão em **Results…**.

## Results…

**Results…** abre a tabela "Clusters between groups": uma linha para cada comparação, par a par quando há mais de dois grupos.

| Coluna | Significado |
|---|---|
| *within each* | O painel, quando a aba tem um |
| Cluster | O nome do cluster |
| Groups | Os grupos comparados, `Control vs Treated` |
| Mean ± SD (n) | A fração média (%) de cada grupo, o desvio-padrão e o número de amostras |
| Test | O teste feito |
| Paired by | A coluna de pareamento, num desenho pareado |
| Statistic | t, F, U, H, W ou χ² |
| p | Desse teste sozinho |
| p adjusted (*método*) | Depois da correção escolhida |
| (sem nome) | `****` < 0,0001, `***` < 0,001, `**` < 0,01, `*` < 0,05, `ns` |
| Note | Por que um teste não pôde ser feito, ou o que saber sobre ele |

Um p em branco é um teste que não pôde ser feito; a Note diz por quê: um grupo só; um grupo com menos de 2 amostras; todas as amostras com o mesmo valor; nenhuma variação dentro dos grupos. Um valor de p nunca é inventado. **Export CSV…** grava `cluster_comparison.csv`; **Copy** põe a tabela na área de transferência.

## Escolher entre as três análises

| Pergunta | Use |
|---|---|
| A fração de um cluster difere, com efeitos claros e amostras suficientes? | **Test** na aba: rápido, aparece no gráfico |
| A fração difere, com clusters raros ou amostras de tamanhos muito diferentes? | Abundância diferencial, [[abundance]]: modela as contagens e achou cerca do dobro das mudanças verdadeiras em simulação |
| As células de um cluster expressam mais de um marcador? | [[differential-expression]] |

Nenhuma delas é um modelo misto: um valor por amostra, comparado por grupo, pareado ou não, opcionalmente ajustado por mais um fator nas duas análises diferenciais.
