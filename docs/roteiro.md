# MIL-C diante do FlowJo e dos outros: o que falta e em que ordem

Levantamento de 22 de novembro de 2026. Serve para decidir o que implementar a seguir, não como
material de divulgação.

## Como ler

- **S** tem · **P** parcial, por plugin ou módulo pago · **N** não tem · **?** não confirmado.
- A coluna **MIL-C** descreve o código de hoje (commit desta data) e foi conferida no próprio
  programa.
- As colunas dos outros programas vêm de uma pesquisa na documentação pública de cada um, com as
  fontes no fim. Onde a página do fabricante não confirmava, ficou **?**. Um **?** não é uma lacuna
  do concorrente: é uma lacuna da pesquisa. Kaluza, FlowLogic e Floreada são os menos documentados.
- "Outros que têm" cita só quem a pesquisa confirmou.

## Comparação por recurso

### Entrada de dados

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| FCS 2.0 / 3.0 / 3.1 | S | S | todos |
| FCS 3.2 | ? | N | FCS Express |
| CSV como eventos | S | N (só exporta) | OMIQ |
| Dados espectrais brutos | S | N | FCS Express, OMIQ |
| Index sort | P (plugin) | N | — |
| Layout de placa | S | N | — |
| Ver e editar keywords | S | N | FCS Express, OMIQ, CellEngine |
| Anotação de amostras em colunas livres | S | **S** | OMIQ, CellEngine |
| Agrupar por metadados, mais de uma coluna | S | **S** | OMIQ, CellEngine |
| Ler workspace do FlowJo (.wsp) | S | **S** (retângulo, polígono, quadrante, booleano e elipse; os três últimos por regras não conferidas com o FlowJo) | FlowKit; Floreada escreve |
| GatingML | N? | N | FCS Express, Cytobank?, FlowKit |

### Compensação e escalas

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| Aplicar a matriz do arquivo ou do workspace | S | **S** | todos |
| Editor de matriz | S | **S** (uma matriz para todas as amostras) | FCS Express, OMIQ, CellEngine, Kaluza |
| Autocompensação por single stains | S | N | FCS Express, OMIQ, FlowLogic, Floreada |
| AutoSpill | S | N | Floreada, R |
| Unmixing espectral | S | N | FCS Express, OMIQ |
| Biexponencial | S | **S** (a do FlowJo, conferida) | FCS Express |
| Logicle, hyperlog | S | N | Floreada (logicle), R |
| Arcsinh, logarítmica | S | **S** | OMIQ, Cytobank, CellEngine |
| Ajuste de escala por parâmetro | S | **S** | quase todos |

### Controle de qualidade

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| Gate manual no tempo | S | P (retângulo sobre Time) | todos |
| flowAI / flowCut / PeacoQC | P (plugins, exigem R) | **S** (método próprio, embutido: vazão e sinal no tempo; não comparado com eles) | OMIQ, FCS Express, Cytobank |
| Relatório de QC | P | P (tabela por amostra do que a limpeza tira) | FCS Express, Kaluza, FlowLogic |

### Gates

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| Retângulo, polígono, faixa (1D) | S | **S** | todos |
| Elipse | S | **S** (desenhada aqui e lida do .wsp) | FCS Express, CellEngine, Kaluza, Floreada |
| Quadrante | S | **S** (idem) | idem |
| Quadrante deslocado, angulado, curvo | S | N | Kaluza (parte) |
| Booleanos (E / OU / NÃO) | S | **S** | FCS Express, Cytobank, Kaluza |
| Autogate por contorno, gate magnético | S | N | Kaluza (contorno) |
| Autogating por modelo | P (plugins) | N | OMIQ, Cytobank, FlowLogic, openCyto |
| Gate sobre clusters | P (plugin) | **S** (cluster do FlowSOM vira população dos arquivos inteiros) | FCS Express, OMIQ, Cytobank |
| Backgating | S | **S** | FCS Express, OMIQ |
| Gate do grupo para todas as amostras | S | **S** | todos |
| Ajuste por amostra | S | **S** | FCS Express, CellEngine, Floreada |

### Gráficos

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| Pontos em pseudocor | S | **S** | todos |
| Contorno | S | **S** | FCS Express, Kaluza |
| Zebra | S | N | — |
| Histograma | S | **S** | todos |
| Sobreposição de amostras (overlay) | S | **S** (histogramas, também escalonados; contornos) | todos |
| Heatmap de populações | S | **S** (clusters) | FCS Express, OMIQ, Cytobank, CellEngine |
| Barras por grupo com cada amostra como ponto | P (table plots) | **S**, com painéis por segunda coluna | FCS Express, OMIQ, Cytobank |
| Violino / box | P (plugin) | **S** | FCS Express, Cytoflow |
| Radar, árvore de fenótipos | N | N | Kaluza |
| Matriz N×N de gráficos | P | **S** (até 12 parâmetros) | FCS Express, Kaluza |
| Cinética no tempo | S | N | FCS Express, Floreada |
| Mapa dividido em painéis por grupo | P (layout) | **S**, por uma ou duas colunas | OMIQ |

### Estatísticas

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| Contagem e % do pai | S | **S** | todos |
| % do avô e do total | S | **S** | todos |
| Mediana, média, média geométrica por população | S | **S** | todos |
| CV, DP, CV robusto, MAD | S | **S** | FCS Express (parte), Kaluza (parte) |
| Percentis arbitrários, moda | S | N | FCS Express (percentil) |
| Stain index | P (plugin) | N | Cytobank |
| Contagem absoluta | P | N | FCS Express |
| Parâmetros derivados e fórmulas | S | N | FCS Express, Kaluza |
| Tabela de estatísticas exportável | S | **S** (CSV e área de transferência) | todos |

### Módulos especializados

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| Ciclo celular | S | N | FCS Express, Kaluza, Floreada |
| Proliferação | S | N | FCS Express |
| Cinética / fluxo de cálcio | S | N | FCS Express, Floreada |
| Comparação de populações (Overton, K-S, χ²) | S | N | — |
| Ensaios com beads | P (plugin) | N | FlowLogic |

### Análise de alta dimensão

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| Downsample | P (plugin) | **S** | FCS Express, OMIQ, Cytobank, CellEngine |
| Plano de amostragem balanceada, à vista antes de concatenar | N | **S** | — (não encontrado em nenhum) |
| Concatenação | S (arquivo novo) | **S** | OMIQ e Cytobank (virtual) |
| UMAP | P (plugin, exige R) | **S**, embutido | FCS Express, OMIQ, Cytobank, CellEngine |
| tSNE | S | **S** | idem |
| opt-SNE / FIt-SNE | ? | P (openTSNE usa FIt-SNE) | OMIQ, Cytobank |
| PCA | P (plugin) | **S** | FCS Express, OMIQ, CellEngine |
| TriMap, PaCMAP, PHATE | P (plugins) | N | OMIQ |
| FlowSOM | P (plugin, exige R) | **S**, embutido | FCS Express, OMIQ, Cytobank, CellEngine |
| PhenoGraph / Louvain / Leiden | P (plugin) | **S** (Louvain, no estilo do PhenoGraph; sem Leiden) | FCS Express, OMIQ, CellEngine |
| X-Shift, SPADE, CITRUS | P / N | N | OMIQ, Cytobank, FCS Express (SPADE) |
| Trajetória (Wanderlust, Wishbone) | P (plugins) | N | OMIQ |
| Explorador de clusters com painéis ligados | S | **S** | FCS Express, OMIQ, FlowLogic |
| Renomear, juntar e anotar clusters | S | **S** | OMIQ, FCS Express |
| Nome automático de clusters | P (plugins) | **S** (sugerido pelos marcadores positivos) | FlowLogic |
| Abundância diferencial (diffcyt, edgeR, SAM) | P (plugin) | **S** (binomial negativa no estilo do edgeR, não conferida com ele) | OMIQ; Cytobank (CITRUS) |
| Normalização de lote (CytoNorm, cyCombine) | P (plugins) | **S** (no estilo do CytoNorm; também em Gates e Statistics; lote alvo) | OMIQ |
| Semente e parâmetros gravados com o resultado | P | **S** | OMIQ (registro de passos) |

### Comparação entre grupos

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| Média e desvio por grupo | P | **S** | vários |
| Teste t (Welch), Mann-Whitney, ANOVA, Kruskal-Wallis | N | **S** | Cytobank, FCS Express (ANOVA), FlowLogic? |
| Correção para múltiplos testes | N | **S** (Benjamini-Hochberg, Holm, Bonferroni) | Cytobank, OMIQ |
| Desenhos pareados ou de medidas repetidas | N | **S** (t pareado, Wilcoxon, ANOVA de medidas repetidas, Friedman) | — |
| Exportar para o Prism | N | **S** (.pzfx, tabela agrupada com réplicas) | FCS Express, OMIQ |

### Figuras, relatórios e exportação

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| Editor de figura por gráfico: tamanho em mm, fonte, cores, legenda | P (limitado; a queixa original) | **S** | — |
| SVG / PDF com texto editável, PNG / TIFF com dpi | P | **S** | OMIQ e Cytobank (SVG, PDF) |
| Estilos de figura salvos e reaplicáveis | N? | **S** | — |
| Rampa de cor à escolha, inclusive cores próprias | P | **S** | — |
| Editor de layout com vários painéis | S | **S** (prancha com letras, largura de coluna de revista) | todos os comerciais |
| Relatório em lote | S | P (PDF: prancha e página de métodos) | todos os comerciais |
| PowerPoint | S | N | FCS Express, OMIQ, FlowLogic |
| Exportar população gatada como FCS / CSV | S | **S** | FCS Express, Cytobank, CellEngine, Floreada |

### Fluxo de trabalho e plataforma

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| Linha de comando para tudo o que a janela faz | P | **S** | R, Python |
| Modelos de análise reaplicáveis | S | P (o arquivo do projeto, por pasta) | todos |
| Pipeline visual de passos | N | N | FCS Express, OMIQ |
| Trilha de auditoria, 21 CFR Part 11 | N | N | Cytobank, CellEngine; FCS Express e OMIQ pagos |
| Colaboração em nuvem | N | N | OMIQ, Cytobank, CellEngine |
| macOS em Apple Silicon, nativo | N (pendente, segundo o fabricante) | **S** | os de navegador |
| Windows | S | **S** | quase todos |
| Linux | ? | P (pelo código, sem instalador) | os de navegador, R, Python |
| Sem internet, sem licença que expira | P (24 h; token de 14 dias) | **S** | FCS Express, Kaluza |
| Gratuito e de código aberto | N | **S** | Cytoflow, FlowKit, R; Floreada é gratuito |

## Onde o MIL-C já está à frente

1. **Figuras.** Tamanho em milímetros, fonte, cores e rótulos por categoria, estilos salvos, texto
   editável no SVG e no PDF. Era a queixa de partida e nenhum dos pesquisados faz mais.
2. **Amostragem balanceada à vista.** O plano mostra o que cada amostra tem, o que vai ceder e
   qual limita, antes de concatenar. Não achei equivalente.
3. **UMAP e FlowSOM sem R.** No FlowJo são plugins que pedem uma instalação de R funcionando.
4. **Subdivisão por duas colunas** no mapa e no gráfico de frequência, sem montar layout à mão.
5. **Registro.** Semente, parâmetros e plano de amostragem ficam com o resultado.
6. **Nativo em Apple Silicon, gratuito, sem licença.**

## Próximas implementações, em ordem

A ordem pesa três coisas: o que um usuário de FlowJo sente falta na primeira hora, o que os
concorrentes fazem melhor que o FlowJo, e o que este experimento em particular precisa.

### Primeira leva: o básico que ainda falta

Sem isto, o programa não substitui o FlowJo no dia a dia, só no fim da análise.

1. ~~**Tabela de estatísticas.**~~ **Feito** (aba Statistics e `milc stats`). Por população e por amostra: contagem, % do pai, do avô e do total,
   mediana, média, média geométrica, CV e CV robusto de qualquer marcador. Exportável, com as
   colunas de anotação junto. É o que alimenta todo gráfico de comparação.
2. ~~**Sobreposição de amostras.**~~ **Feito** (aba Gates, com várias amostras selecionadas). Histogramas e contornos de várias amostras ou grupos num mesmo
   gráfico, com deslocamento. É como se confere um FMO ou um controle.
3. ~~**Gates que faltam.**~~ **Feito**: elipse, quadrante e booleanos (E / OU / NÃO), desenhados no
   programa. Ler elipses e quadrantes de um workspace do FlowJo continua por fazer.
4. ~~**Escala por parâmetro.**~~ **Feito**, menos a logicle (biexponencial ajustável, arcsinh,
   logarítmica e linear). Ajustar a biexponencial (largura, décadas negativas), e oferecer
   logicle e arcsinh também na visualização.
5. ~~**Ver e editar a matriz de compensação.**~~ **Feito** (uma matriz para todas as amostras).
6. ~~**Exportar uma população gatada**~~ **Feito**: FCS ou CSV, por amostra.
7. ~~**Gráfico de contorno e densidade.**~~ **Feito** (contornos de probabilidade; zebra não).

### Segunda leva: onde vale passar o FlowJo

São as lacunas do FlowJo que mais de um concorrente já cobre, e todas servem ao objetivo do MIL-C.

8. ~~**Testes entre grupos.**~~ **Feito** (aba Clusters e aba Statistics). t, Mann-Whitney, ANOVA e Kruskal-Wallis sobre a frequência dos
   clusters e sobre as estatísticas da tabela, com correção para múltiplos testes e o valor de p
   no gráfico. O FlowJo não tem.
9. ~~**Saída para o Prism.**~~ **Feito** (arquivo `.pzfx`; falta conferir a abertura no Prism). Tabela já no formato de colunas agrupadas.
10. **Trabalhar os clusters.** Renomear, juntar, anotar; transformar um cluster ou um laço em
    população; sugerir nome pelo perfil de marcadores.
11. **Mais um método de cluster e um de redução.** PhenoGraph / Leiden, que não pede o número de
    clusters, e PCA.
12. ~~**Normalização de lote.**~~ **Feito**, no estilo do CytoNorm, sem comparação com o pacote R.
    Este experimento foi adquirido em três dias (9, 13 e 19 de novembro),
    e parte da diferença entre dias pode ser do aparelho. CytoNorm é o método usual.
13. **Limpeza automática.** PeacoQC ou flowAI, e o gráfico de sinal no tempo.
14. **Abundância diferencial.** diffcyt / edgeR, com volcano.
15. **Prancha de vários painéis** e relatório em lote, com os estilos de figura que já existem.

### Terceira leva: especializações

16. Violino e box por marcador; matriz N×N.
17. Backgating.
18. Proliferação, ciclo celular e cinética.
19. Autocompensação por single stains e AutoSpill.
20. GatingML, e escrever workspace do FlowJo.
21. Autogating por modelo.
22. Unmixing espectral.
23. Instalador para Linux; trilha de auditoria.

### O que eu não faria

- **Colaboração em nuvem.** É o produto de OMIQ, Cytobank e CellEngine, e não combina com um
  programa de mesa, gratuito e que roda sem internet.
- **CITRUS, SPADE, X-Shift.** Métodos que a área vem trocando por FlowSOM e Leiden.

## Onde estamos (23 de novembro de 2026)

Feito desde o levantamento: os itens 1 a 9 e o 12 — tabela de estatísticas, sobreposição de
amostras, gates de elipse, quadrante e booleanos, escala por parâmetro, edição da compensação,
exportação de população, contorno, testes entre grupos, saída para o Prism e normalização de
lote.

Feito em seguida, da lista "Próximas implementações" que estava aqui: trabalhar os clusters
(nomes, notas, junção, nome sugerido, cluster como população); prancha de vários painéis e
relatório em PDF; limpeza automática da aquisição com gráfico no tempo; PhenoGraph e PCA;
abundância diferencial com volcano; testes pareados, de medidas repetidas e ANOVA de Welch;
normalização também em Gates e Statistics, com lote alvo; leitura de quadrantes, booleanos e
elipses do workspace do FlowJo; violino e box, matriz N×N e backgating.

### A conferir, antes de qualquer item novo

Coisas feitas que ainda dependem de um olhar de fora:

- **Abrir um `.pzfx` no Prism.** É o único jeito de saber se o Prism aceita o arquivo.
- **Instalar o `.msi` num Windows de verdade.** O executável foi testado; o instalador, não.
- **Usar a janela com mouse e teclado.** Desenho de gates, laço, seleção múltipla e os diálogos
  foram testados por código e por capturas, não à mão.
- **Uma estatística lado a lado com o FlowJo.** As definições são as dele; os números não foram
  comparados.
- **Um workspace com elipse, quadrante ou booleano.** As regras de leitura vieram do FlowKit e
  do CytoML; nenhum arquivo nosso tem esses gates. A contagem do FlowJo aparece ao lado de cada
  população lida assim: é só comparar.
- **A abundância diferencial ao lado do edgeR/diffcyt**, e o PhenoGraph ao lado do pacote. Os
  dois foram calibrados por simulação, não comparados com o original.
- **A limpeza automática numa amostra que se sabe ruim.** Em dados simulados acerta o trecho; nos
  arquivos reais tirou de 0 a 8%, sem gabarito. A `14Ctrl1` perde os primeiros ~23 s.
- **A normalização "em Gates e Statistics".** Os gates foram desenhados nos valores originais:
  depois de ligá-la, vale rever cada um.

### Próximas implementações, na ordem em que eu faria

1. **Comparar com os originais** o que foi escrito a partir dos artigos: abundância diferencial
   com o edgeR, PhenoGraph com o pacote, limpeza com o PeacoQC, normalização com o CytoNorm
   inteiro. Precisa de R com Bioconductor.
2. **Logicle**, e matriz de compensação por amostra.
3. **Expressão diferencial por cluster** (o "DS" do diffcyt): a mediana de um marcador dentro de
   cada cluster, entre grupos.
4. **Leiden** no lugar do Louvain, e tamanho mínimo de cluster no PhenoGraph, que hoje divide
   demais as populações grandes na resolução 1.
5. **Prancha guardada com o projeto**, e mais de uma prancha; hoje ela dura a sessão.
6. **Limpeza pela linha de comando**, e o resumo dela no relatório.
7. **Manual embutido**, em inglês e português, como no MIL-Q e no MIL-X.
8. **Distribuição.** Publicar a versão 0.1.0; instalador para Linux; assinatura dos instaladores
   (hoje o macOS e o Windows avisam na primeira abertura).

Mais adiante, sem mudança: proliferação, ciclo celular e cinética; autocompensação e AutoSpill;
GatingML; autogating; unmixing espectral; trilha de auditoria.

## Fontes

Documentação e páginas dos fabricantes, consultadas em 22 de novembro de 2026:

- **FlowJo:** docs.flowjo.com (graphs-and-gating, advanced-gates, gw-transform-overview,
  workspaces-and-samples, ws-statdefinitions, experiment-based-platforms, advanced-features,
  graphical-reports, tabular-reports, plugins-2) e flowjo.com/exchange.
- **FCS Express:** manual em fcsexpressdownloads.s3.amazonaws.com (new_features),
  denovosoftware.com (features/gates, system-requirements, full-access).
- **OMIQ:** omiq.ai/features.
- **Cytobank:** notas de versão em beckman.com/flow-cytometry/software/cytobank-premium.
- **CellEngine:** docs.cellengine.com (algorithms, scales, pricing).
- **Kaluza:** beckman.com/flow-cytometry/software/kaluza (features, learning-center).
- **FlowLogic:** inivai.com/flowlogic.
- **Infinicyt:** bdbiosciences.com (infinicyt-software).
- **Floreada:** floreada.io.
- **Código aberto:** cytoflow.readthedocs.io, flowkit.readthedocs.io.
- **Revisão cruzada:** graphpad.com/guides/best-flow-cytometry-analysis-software — publicada pelo
  fabricante do Prism, do OMIQ e do FCS Express, portanto não neutra.

Não conferidos nesta pesquisa: Flowing Software, Cytek SpectroFlo, Tercen, WinList / ModFit,
NovoExpress, CytExpert, e os pacotes de R (flowCore, openCyto, CATALYST, Spectre, CytoExploreR),
citados de conhecimento prévio.
