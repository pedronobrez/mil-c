---
title: A aba Statistics
---

A aba Statistics é uma tabela: uma linha para cada população de cada amostra, com as contagens e porcentagens que toda linha tem e, para os parâmetros pedidos, as estatísticas que um citometrista relata. As colunas das próprias amostras, vindas da [[samples-tab|aba Samples]] (Group, dia, doador…), ficam ao lado do nome da amostra, de modo que a tabela pode ser ordenada aqui e comparada por grupo aqui ou em qualquer outro lugar.

Ela lê os arquivos individuais com os seus gates, não os eventos concatenados: todo evento da população é contado, não só os sorteados para a [[concatenation|concatenação]].

## Escolher o que calcular

O painel da esquerda tem três listas e os botões.

- **Populations** — a árvore de gates, com uma caixa de marcação em cada população. Só aparecem as populações que toda amostra com gates tem, porque cada uma precisa de uma linha em todas as amostras. Um gate desenhado desde a última visita entra marcado; um que você desmarcou continua desmarcado. "All events" é a raiz, o arquivo sem gate.
- **Parameters** — todos os parâmetros das amostras, pelo nome do marcador onde há um. A fluorescência é tomada compensada e nas unidades do próprio arquivo, como o FlowJo relata: a coluna lida é a `Comp-` sempre que o parâmetro é compensado (ver [[compensation]]). Com o [[spectral-unmixing|unmixing espectral]] aplicado, aparecem os fluorocromos desmisturados no lugar dos detectores que eles substituem.
- **Statistics** — marque qualquer uma das oito abaixo. O que você marca é lembrado na próxima vez; na primeira, vem marcada a Median.
- **Only the samples ticked in Samples** — marcado por padrão. Desmarque para ter uma linha para cada amostra utilizável.
- **Compute** — monta a tabela, em segundo plano; a barra de status mostra qual amostra está sendo lida.

A contagem e as três porcentagens são sempre calculadas, mesmo sem nenhum parâmetro marcado.

## As estatísticas

Toda estatística é calculada sobre os valores dos eventos como o arquivo os guarda, depois da compensação e antes de qualquer escala de exibição: uma mediana de fluorescência é uma fluorescência, não uma posição num eixo biexponencial. As definições são as do FlowJo, de modo que um número daqui pode ser posto ao lado de um de lá.

| Estatística | Definição |
|---|---|
| Median | O valor do meio dos eventos da população |
| Mean | A média aritmética |
| Geometric mean | exp(média(ln x)) sobre os eventos com x > 0. Eventos em zero ou abaixo não têm logaritmo e ficam de fora |
| SD | Desvio-padrão, com n − 1 |
| CV | 100 × SD / Mean |
| Robust SD | Metade da distância entre os percentis 15,87 e 84,13: o desvio-padrão de uma normal com essa largura, surdo a valores extremos |
| Robust CV | 100 × Robust SD / Median |
| MAD | Mediana dos desvios absolutos em relação à mediana |

Uma população sem eventos numa amostra dá uma célula vazia, não um zero. Um CV cuja média (ou mediana) é zero também fica vazio.

## Colunas

| Coluna | Significado |
|---|---|
| Sample | O nome da amostra |
| Suas colunas | Uma para cada coluna da tabela de Samples, como digitada lá |
| Population | O nome da população; "All events" para a raiz |
| Path | O caminho completo do gate. Só aparece quando duas populações têm o mesmo nome |
| Count | Eventos na população |
| % of parent | Contagem sobre a da população mãe, × 100 |
| % of grandparent | Contagem sobre a mãe da mãe, × 100 (vazia no primeiro nível de gates) |
| % of total | Contagem sobre todos os eventos do arquivo, × 100 |
| *Marcador* · *Estatística* | Uma coluna para cada parâmetro e estatística marcados, por exemplo `CD4 · Median` |

Clique num cabeçalho para ordenar por ele; números são ordenados como números. Contagens aparecem com separador de milhar, porcentagens e CVs com duas casas decimais, os demais números com três algarismos significativos (uma casa decimal a partir de 100). Os valores exportados e copiados guardam seis algarismos significativos.

A legenda acima da tabela diz quantas linhas, amostras e populações ela tem, e se as amostras foram normalizadas por lote ao serem lidas (ver [[normalisation]]). Quando os gates, as amostras ou suas colunas mudam depois do Compute, aparece um aviso ao lado: **Computed before the last change — Compute again**. A tabela continua mostrando o que foi calculado até você calcular de novo.

## Levar a tabela para outro lugar

- **Export CSV…** — a tabela inteira como `statistics.csv` (ou o nome que você escolher). O arquivo é UTF-8 com marca de ordem de bytes, para que o Excel no Windows leia `CD4 · Median` corretamente.
- **Copy** — a tabela inteira, separada por tabulações e com o cabeçalho, para colar no Excel ou no Prism.
- **Prism…** — uma coluna da tabela como arquivo do GraphPad Prism: ver [[tables-and-prism]].
- **Compare groups…** — testa uma coluna entre grupos de amostras: abaixo.

A mesma tabela pode ser feita sem a janela: `milc stats` na [[command-line|linha de comando]].

## Compare groups

**Compare groups…** abre um diálogo que pergunta o que comparar e como.

| Campo | O que faz |
|---|---|
| **Compare** | A coluna da tabela a testar: Count, qualquer porcentagem ou qualquer coluna de estatística. Começa em % of parent |
| **Between** | As colunas das amostras que formam os grupos. Marque uma (Group) ou várias: então cada combinação é um grupo, escrito `Treated · 7` |
| **Within each** | Opcional. Uma comparação separada para cada valor desta coluna — os grupos dentro de cada dia, por exemplo |
| **Paired by** | Opcional. A coluna que diz de qual animal ou doador é cada amostra, quando os mesmos foram medidos em todos os grupos |
| **Test** | **Welch t test / ANOVA**, **Mann-Whitney / Kruskal-Wallis** ou **Welch t test / Welch ANOVA** |
| **Adjust p by** | **Benjamini-Hochberg**, **Holm**, **Bonferroni** ou **None** |
| **With more than two groups, every pair as well** | Marcado por padrão: além do teste de todos os grupos juntos, um teste para cada par |

É feito um teste para cada população (e cada valor de *Within each*), e os valores de p são ajustados sobre todos os testes juntos, pares incluídos. A unidade é a amostra, nunca a célula: cinco animais são cinco observações, quantas células cada um tenha dado. Os testes, o pareamento e as notas estão explicados em [[cluster-statistics]], que usa a mesma maquinaria para os clusters.

O resultado abre numa janela própria, uma tabela com o nome do valor comparado:

| Coluna | Significado |
|---|---|
| *Within each* | O valor dessa coluna, quando uma foi escolhida |
| Population | O que foi testado |
| Groups | Os grupos deste teste, `Control vs Treated` |
| Mean ± SD (n) | A média, o desvio-padrão e o número de amostras de cada grupo |
| Test | O teste feito: Welch t, ANOVA, Welch ANOVA, Mann-Whitney, Kruskal-Wallis, Paired t, Wilcoxon signed-rank, RM ANOVA (Greenhouse-Geisser) ou Friedman |
| Paired by | A coluna de pareamento, num desenho pareado |
| Statistic | t, F, U, H, W ou χ² do teste |
| p | O p desse teste sozinho |
| p adjusted (*método*) | Depois do ajuste pelo número de testes |
| (sem nome) | Estrelas pelo p ajustado: `****` < 0,0001, `***` < 0,001, `**` < 0,01, `*` < 0,05, `ns` no resto |
| Note | Por que um teste não pôde ser feito, ou o que saber sobre ele |

Um p em branco é um teste que não pôde ser feito; a Note diz por quê. Valores de p abaixo de 0,0001 aparecem como `< 0.0001`. **Export CSV…** (como `group_comparison.csv`), **Copy** e **Close** ficam abaixo da tabela.

## Limites

- As estatísticas seguem as definições do FlowJo, mas não foram postas ao lado dos números do FlowJo nos mesmos arquivos (ver [[measured-facts]]).
- Os testes usam um valor por amostra, por grupo, pareados ou não. Não há modelos mistos.
- Uma população que alguma amostra com gates não tem nem é oferecida; dê o gate a essa amostra primeiro.
