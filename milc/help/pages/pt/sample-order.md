---
title: A ordem dos grupos
---

Qual grupo vem primeiro é um detalhe que aparece em toda parte: as barras de um gráfico de
frequência, as entradas de uma legenda, os painéis de um mapa dividido, as linhas de uma tabela,
as colunas de um arquivo do Prism. O MIL-C guarda uma ordem para cada coluna da tabela da aba
Samples e a usa em todos esses lugares, para que Control fique à esquerda de Treated no mapa, no
gráfico e na tabela exportada.

## Ordem de contagem

Sem intervenção, os valores vêm na ordem em que uma pessoa conta, não na ordem em que um
computador ordena texto:

| Valores | Ordem de contagem | Ordem de texto simples |
|---|---|---|
| dia 3, 7, 14 | 3, 7, 14 | 14, 3, 7 |
| `sample_2`, `sample_10` | `sample_2`, `sample_10` | `sample_10`, `sample_2` |
| Treated, control | control, Treated: maiúsculas não importam | Treated, control: maiúsculas primeiro |

Números dentro de um valor são comparados como números, e letras maiúsculas e minúsculas como a
mesma letra. Os arquivos de uma pasta são listados na mesma ordem.

## Dando uma ordem

**Order…** abre o diálogo **Order of the groups**. Ele está na aba Samples, acima da tabela, e na
aba Clusters, onde costuma ser mais procurado.

- **Column** escolhe a coluna cujos valores aparecem. O diálogo abre na primeira coluna, ou naquela
  para a qual o botão foi apertado.
- Arraste um valor para cima ou para baixo na lista, ou selecione-o e use **Move up** e **Move
  down**.
- **Reverse** inverte a lista.
- **Counting order** volta à ordem que ninguém escolheu: 3, 7, 14; A, B, C.
- A linha abaixo da lista diz qual é a situação: `In the order given here.` (na ordem dada aqui),
  `In counting order: nobody chose another.` (em ordem de contagem) ou `No values in this column
  yet.` (a coluna ainda não tem valores).

Cada coluna guarda a sua ordem. Dá para arrumar uma coluna, trocar para outra e arrumar essa, e
tudo o que foi arrumado é aplicado junto no **OK**. **Cancel** não muda nada.

## Onde vale

Depois de aplicada, a ordem é usada por toda parte do programa que mostra grupos:

- os gráficos de frequência, violinos e caixas, e as legendas de todo gráfico;
- os painéis de um mapa dividido por uma coluna, em [[explore-tab]] e [[clusters-tab]];
- os grupos dos testes entre grupos ([[statistics-tab]], [[cluster-statistics]]);
- as tabelas de um arquivo do Prism ([[tables-and-prism]]);
- as curvas do diálogo de normalização de lotes ([[normalisation]]) e as opções de
  **Reference** e **Bring to**;
- figuras e pranchas feitas a partir de qualquer um desses ([[figures]]).

Nada é recalculado quando a ordem muda: a ordem é do que se mostra, não do que se mediu. Todo
gráfico que mostra grupos é simplesmente redesenhado, e a barra de status diz `The order of the
groups changed`.

## Grupos feitos de duas colunas

Quando **Group by** marca duas colunas, um grupo é uma combinação, escrita com um ponto:
`Treated · 7`. Um rótulo assim é ordenado pela primeira coluna e depois pela segunda, cada uma na
sua ordem. Com Group na ordem Control, Treated e Day em ordem de contagem, os grupos vêm como
`Control · 3`, `Control · 7`, `Control · 14`, `Treated · 3`… Qual coluna fica por fora segue a
ordem das colunas na tabela da aba Samples, qualquer que seja a ordem em que foram marcadas.

## Regras que vale saber

- Um valor tem um só lugar, em qualquer coluna em que esteja. Se duas colunas tiverem o mesmo
  valor e lhe derem lugares diferentes, vale a ordem da primeira coluna.
- Valores não citados na ordem de uma coluna — um grupo novo digitado depois na tabela — vêm
  depois dos citados, em ordem de contagem.
- Renomear uma coluna mantém a sua ordem; remover uma coluna remove a sua ordem.
- A ordem fica no arquivo de projeto da pasta ([[projects-and-saving]]) e está lá da próxima vez
  que a pasta for aberta.
