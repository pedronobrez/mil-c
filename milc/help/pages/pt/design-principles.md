---
title: Princípios de projeto
---

O MIL-C foi escrito para substituir um caminho pelo FlowJo — fazer gates, reduzir a amostra, concatenar, rodar os plugins de UMAP, FlowSOM e ClusterExplorer — e para fazer três coisas que esse caminho fazia mal: figuras que podem ser publicadas como estão, amostras que pesam o mesmo na análise e um registro de como cada resultado foi feito. As escolhas abaixo decorrem disso. Elas explicam por que o programa se comporta como se comporta onde difere do que você talvez espere.

## O balanceamento é decidido antes de qualquer sorteio

Uma amostra tratada muitas vezes traz várias vezes as células do seu controle, e um mapa feito com tudo é um mapa do grupo tratado. Por isso o plano de amostragem fica à vista antes de qualquer concatenação: o que cada amostra tem, o que vai ceder e qual amostra definiu o limite. Há três modos de balancear — o mesmo número de cada amostra, o mesmo total de cada grupo, ou um teto por amostra — e um mínimo abaixo do qual uma amostra fica de fora em vez de puxar todas as outras para baixo. O custo do balanceamento é lido antes de ser pago. Ver [[concatenation]].

## Os gates são a única definição de uma população

Uma população é o que os seus gates selecionam, e toda parte do programa pergunta à árvore de gates: as [[statistics-tab|estatísticas]], a exportação de eventos, a concatenação, a limpeza do [[quality-control|controle de qualidade]] (que vira uma primeira população, Clean, de onde todo o resto parte) e até um cluster, que vira população por um gate que situa cada evento no mapa treinado do FlowSOM (ver [[clusters-tab]]). Um gate desenhado numa amostra é de todas, a menos que você diga o contrário. Os gates vêm do workspace do FlowJo ou são desenhados aqui, ficam guardados ao lado dos arquivos e saem de novo como workspace ou Gating-ML (ver [[workspaces-and-gating-ml]]).

## Nada muda pelas suas costas

- Os arquivos `.fcs` são só lidos. Exportações, arquivos concatenados e projetos são gravados como arquivos novos.
- O resultado de um algoritmo nunca é sobrescrito. Nomes, notas, cores e junções de clusters ficam por cima dele, de modo que uma junção pode ser desfeita e nada do que o algoritmo disse se perde.
- Quando os valores por baixo de algo mudam depois — outra matriz de compensação, amostras normalizadas por lote, um gate movido — o programa avisa, em vez de recalcular em silêncio ou manter o antigo em silêncio: a aba Gates marca os gates desenhados sobre os valores antigos até você apertar **Reviewed**, uma linha acima das abas diz que os eventos concatenados estão desatualizados, e a tabela de estatísticas diz **Computed before the last change — Compute again**.
- Um gate proposto pelo programa é ajustado de novo quando a sua população muda; um gate desenhado por você, ou um proposto que você depois moveu, não é sobrescrito — uma barra oferece **Refit** ou **Keep**.
- Um valor de p nunca é inventado. Onde um teste não pode ser feito, a célula fica em branco e a nota diz por quê.
- A ordem dos grupos é uma ordem de exibição. Mudá-la redesenha; nada é recalculado (ver [[sample-order]]).

## Números nas unidades em que foram medidos

As estatísticas são de valores compensados nas unidades do próprio arquivo, pelas definições do FlowJo: uma mediana de fluorescência é uma fluorescência, não uma posição num eixo de exibição. As escalas de exibição são escolhidas onde o dado é mostrado ou entregue a um algoritmo, e os eventos concatenados são guardados sem transformação. Ver [[scales-and-transforms]].

## A unidade é a amostra

Os testes entre grupos tomam um valor por amostra. Cinco animais são cinco observações, quantas células cada um tenha dado; testar células contra células acharia toda diferença significativa e não diria nada sobre os animais. Onde muitas coisas são testadas ao mesmo tempo — cada cluster, cada população — os valores de p são ajustados por isso. Ver [[cluster-statistics]].

## Os padrões do FlowJo, para que um primeiro resultado seja comparável

Os algoritmos recebem o que o FlowJo entrega aos seus plugins — marcadores na escala biexponencial, de 0 a 1023 — e rodam com os padrões dos plugins: UMAP com 15 vizinhos e distância mínima 0,5, tSNE com perplexidade 30, FlowSOM numa grade 10 × 10. Um primeiro mapa pode então ser posto ao lado do que ele substitui. Ver [[embeddings]] e [[clustering-methods]].

## Medido, não alegado

Os métodos foram escritos a partir dos artigos e depois postos ao lado dos programas que os tornaram conhecidos, com a mesma entrada. O que foi medido está escrito, inclusive onde o MIL-C difere e onde nada foi conferido; a comparação achou defeitos duas vezes, e eles foram corrigidos. Ver [[measured-facts]].

## Um registro

Sementes, configurações, marcadores, o plano de amostragem, a normalização e os limites do controle de qualidade são salvos com o resultado, e o relatório das pranchas os declara numa página própria. A mesma semente e as mesmas configurações dão o mesmo resultado. Ver [[projects-and-saving]].

## O núcleo não precisa da janela

A análise — ler arquivos, gates, amostragem, os algoritmos, as estatísticas, as figuras — é um pacote Python que nunca importa o Qt. A janela é um modo de conduzi-lo; a [[command-line|linha de comando]] é outro, e faz tudo o que a janela faz sem ela. As mesmas funções podem ser chamadas de um script Python. Isso mantém os métodos testáveis por si mesmos e os resultados reprodutíveis sem cliques.

## Sem R, sem licença, sem rede

UMAP, FlowSOM, PhenoGraph, abundância e expressão diferenciais, controle de qualidade e normalização estão escritos em Python, então não é preciso instalar o R, ao contrário dos plugins do FlowJo. O programa é gratuito, funciona sem internet e é nativo no Apple Silicon. Colaboração em nuvem ficou de fora de propósito: não combina com um programa de mesa que roda sem internet.

## As figuras são o produto

Todo gráfico vira uma figura em milímetros, na fonte e nas cores que uma revista pede, com texto editável no SVG e no PDF, e a prévia é desenhada pelo mesmo código que grava o arquivo. Vários painéis vão numa página só em [[the-plate]]. Ver [[figures]].

## Uma língua na interface, um manual traduzido

A interface é em inglês e continua em inglês, para que um controle tenha o mesmo nome em todo laboratório, captura de tela e seção de métodos de artigo. Este manual também está em português; aqui, os nomes dos controles ficam em inglês e são explicados.
