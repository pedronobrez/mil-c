---
title: Controle de qualidade da aquisição
---

Um entupimento, uma bolha de ar, uma mudança de pressão ou um laser que deriva estragam os eventos
adquiridos enquanto duram. O MIL-C procura esses trechos como o flowAI e o PeacoQC procuram — na
ordem em que os eventos foram adquiridos, um número fixo de eventos por vez — e o que ele acha é
uma lista de intervalos de tempo. Aplicados, eles viram uma população, **Clean**, no topo da
árvore de gates de cada amostra: todo gate, estatística e concatenação passa a partir dos eventos
limpos. Não é preciso R.

## Abrindo

**Quality…** fica na [[gates-tab]], entre as ferramentas de gate; escolha uma amostra antes. O
diálogo **Quality of the acquisition: amostra** abre naquela amostra, já verificada com os limites
aplicados da última vez (ou os padrões). Nada é retirado só por olhar.

## O que é verificado

Os eventos são postos na ordem do parâmetro **Time** e cortados em bins com o mesmo número de
eventos. Um bin é suspeito quando:

| Verificação | O que acha |
|---|---|
| **Flow rate** (vazão) | os eventos por segundo do bin estão longe da vazão usual da aquisição: um entupimento, uma bolha. Uma vazão dentro de 1,5 vez a usual nunca é defeito |
| **Signal** (sinal) | a mediana móvel, de cinco bins, das medianas de uma população fica a mais que o limite de MAD da mediana dessa população no arquivo inteiro: um degrau, qualquer coisa abrupta |
| **Drift** (deriva) | o nível lento de uma população, seguido numa janela longa, se afasta de onde fica por mais tempo: um desvio moderado que dura |
| **Shift** (desvio curto) | um trecho curto, de poucos bins, cujo nível se destaca dos bins dos dois lados: um desvio moderado que não dura |

Depois, um trecho de menos de **Consecutive bins** bins bons entre dois trechos removidos também é
removido.

Todo parâmetro menos Time é julgado, a fluorescência compensada e cada um na escala em que o seu
eixo é desenhado. Um parâmetro com várias populações — um marcador que metade das células tem — é
cortado nos vales do seu histograma e seguido uma população por vez, para que um bin com uma
mistura um pouco diferente de células positivas e negativas não seja tomado por defeito. Uma
população pequena demais para dar uma mediana na maioria dos bins (menos de 20 eventos em metade
deles) não é seguida. Um defeito só é achado enquanto atinge menos da metade da aquisição: com
mais, a parte ruim vira a referência.

## Os limites

| Controle | Padrão | Significado |
|---|---|---|
| **MAD limit** | 6 | a quantos MADs do nível da própria aquisição o sinal, ou a vazão, precisa estar para um bin ser suspeito. Maior remove menos |
| **Consecutive bins** | 5 | bins bons entre dois trechos removidos só ficam quando são pelo menos esses |
| **Events per bin** | auto | auto é 1/500 do arquivo e pelo menos 150 eventos. Bins maiores são mais estáveis e veem menos de um defeito curto |
| **Drift limit** | 0,25 | o menor desvio duradouro tomado por deriva, como fração do espalhamento dos eventos do parâmetro. **off** (0) não procura derivas |
| **Shift limit** | 0,25 | o menor afastamento de um trecho curto tomado por desvio, na mesma unidade; ele precisa se destacar também pelo limite de MAD. **off** não procura |

Mude um limite e **Check again** verifica esta amostra com os novos limites — e todas as amostras,
se elas tiverem sido verificadas com os antigos.

## Lendo os gráficos

A linha do alto resume o que sairia: `To remove: 3.2% of 120,000 events, in 2 stretch(es) of time
— flow rate 1,200 events; signal 2,640 (CD4, CD8).` Quando mais de um quinto do arquivo sairia,
ela acrescenta: `look at the plots before applying.` (olhe os gráficos antes de aplicar).

- **Gráfico de cima**: a vazão de cada bin ao longo do tempo (em segundos quando o arquivo traz
  `$TIMESTEP`, senão nas unidades do próprio arquivo). Linhas tracejadas: até onde ela pode ir;
  pontos vermelhos: bins além disso.
- **Gráfico de baixo**, para a população escolhida em **Signal of**: a mediana de cada bin, com
  os limites tracejados e os bins suspeitos em vermelho; uma linha grossa, o nível lento, vermelha
  onde derivou além dos limites pontilhados; barras vermelhas sobre trechos curtos que se
  destacaram dos vizinhos (barras tracejadas).
- **Sombreado** nos dois: os trechos a remover.

**Signal of** começa na população que deu errado por mais tempo. Os eixos de tempo dos dois
gráficos andam juntos.

## Todas as amostras

**Check every sample** verifica cada arquivo com os limites atuais, fora da linha de execução da
janela, e preenche uma tabela: **Sample**, **Events**, **To remove** (%), os eventos achados por
**Flow rate**, **Signal**, **Drift** e **Shift**, e uma **Note** com os parâmetros culpados — ou o
motivo de o arquivo não ter podido ser verificado. Uma amostra verificada tem uma caixa, marcada.

Um arquivo não pode ser verificado quando não tem parâmetro Time, tem eventos de menos (são
precisos 20 bins: 3.000 eventos no menor bin), tem Time faltando em alguns eventos, tem o mesmo
Time para todos os eventos, ou tem Time andando para trás (várias aquisições num arquivo, ou um
relógio zerado).

## Aplicando e removendo

- **Apply to this sample** remove os trechos achados na amostra mostrada.
- **Apply to every sample** (depois de **Check every sample**) remove de cada amostra marcada os
  seus próprios trechos; uma amostra desmarcada fica com todos os eventos.
- **Remove** tira **Clean** da árvore de todas as amostras; os gates voltam para onde estavam.

Aplicar põe **Clean** no topo da árvore de todas as amostras — as que ficam inteiras também a
recebem, com todos os seus eventos, para que as árvores continuem iguais — e todo gate existente
passa para baixo dela. A população escolhida na aba Samples acompanha (`Cells/Live` vira
`Clean/Cells/Live`). Gates automáticos ajustados aos dados são reajustados sobre os eventos
limpos ([[auto-gates]]). A linha no pé do diálogo diz o que as árvores têm e os limites aplicados
por último; eles são salvos com o projeto e escritos no relatório de [[the-plate]].

## Comparado com os originais

Rodado ao lado do PeacoQC 1.22.0 e do flowAI 1.42.0 sobre os mesmos valores (detalhes em
[[measured-facts]]):

- Em defeitos conhecidos, escritos num arquivo limpo, e em defeitos simulados, o MIL-C removeu
  inteiros um degrau e um entupimento; o PeacoQC, que não verifica a vazão, não achou os
  entupimentos.
- Em 32 arquivos reais, o MIL-C removeu 5,25% de todos os eventos, o PeacoQC 5,13% e o flowAI
  11,97%.
- Os limites de deriva e de desvio curto foram medidos em cópias de arquivos reais com desvios de
  tamanho conhecido: os arquivos reais perdem o mesmo com **Drift limit** entre 0,20 e 0,33 e
  **Shift limit** entre 0,20 e 0,30; 0,25 é o meio das duas faixas.

## Pela linha de comando

`milc qc DADOS --out qc.csv --plots qc.pdf` verifica cada arquivo e grava a tabela e uma página de
gráficos por amostra, sem remover nada. `--qc` em `concat`, `stats`, `export` e `gates` remove os
trechos só naquela execução, sem gravar o projeto. Veja [[command-line]].
