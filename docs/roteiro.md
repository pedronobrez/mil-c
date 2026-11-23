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
| Ler workspace do FlowJo (.wsp) | S | **S** (retângulo, polígono, compensação) | FlowKit; Floreada escreve |
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
| flowAI / flowCut / PeacoQC | P (plugins, exigem R) | N | OMIQ, FCS Express, Cytobank |
| Relatório de QC | P | N | FCS Express, Kaluza, FlowLogic |

### Gates

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| Retângulo, polígono, faixa (1D) | S | **S** | todos |
| Elipse | S | **S** (desenhada aqui; não lida do .wsp) | FCS Express, CellEngine, Kaluza, Floreada |
| Quadrante | S | **S** (idem) | idem |
| Quadrante deslocado, angulado, curvo | S | N | Kaluza (parte) |
| Booleanos (E / OU / NÃO) | S | **S** | FCS Express, Cytobank, Kaluza |
| Autogate por contorno, gate magnético | S | N | Kaluza (contorno) |
| Autogating por modelo | P (plugins) | N | OMIQ, Cytobank, FlowLogic, openCyto |
| Gate sobre clusters | P (plugin) | P (laço com estatísticas, sem virar população) | FCS Express, OMIQ, Cytobank |
| Backgating | S | N | FCS Express, OMIQ |
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
| Violino / box | P (plugin) | N | FCS Express, Cytoflow |
| Radar, árvore de fenótipos | N | N | Kaluza |
| Matriz N×N de gráficos | P | N | FCS Express, Kaluza |
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
| PCA | P (plugin) | N | FCS Express, OMIQ, CellEngine |
| TriMap, PaCMAP, PHATE | P (plugins) | N | OMIQ |
| FlowSOM | P (plugin, exige R) | **S**, embutido | FCS Express, OMIQ, Cytobank, CellEngine |
| PhenoGraph / Louvain / Leiden | P (plugin) | N | FCS Express, OMIQ, CellEngine |
| X-Shift, SPADE, CITRUS | P / N | N | OMIQ, Cytobank, FCS Express (SPADE) |
| Trajetória (Wanderlust, Wishbone) | P (plugins) | N | OMIQ |
| Explorador de clusters com painéis ligados | S | **S** | FCS Express, OMIQ, FlowLogic |
| Renomear, juntar e anotar clusters | S | N | OMIQ, FCS Express |
| Nome automático de clusters | P (plugins) | N | FlowLogic |
| Abundância diferencial (diffcyt, edgeR, SAM) | P (plugin) | N | OMIQ; Cytobank (CITRUS) |
| Normalização de lote (CytoNorm, cyCombine) | P (plugins) | N | OMIQ |
| Semente e parâmetros gravados com o resultado | P | **S** | OMIQ (registro de passos) |

### Comparação entre grupos

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| Média e desvio por grupo | P | **S** | vários |
| Teste t, Mann-Whitney, ANOVA | N | N | Cytobank, FCS Express (ANOVA), FlowLogic? |
| Correção para múltiplos testes | N | N | Cytobank, OMIQ |
| Exportar para o Prism | N | N (CSV genérico) | FCS Express, OMIQ |

### Figuras, relatórios e exportação

| Recurso | FlowJo | MIL-C | Outros que têm |
|---|---|---|---|
| Editor de figura por gráfico: tamanho em mm, fonte, cores, legenda | P (limitado; a queixa original) | **S** | — |
| SVG / PDF com texto editável, PNG / TIFF com dpi | P | **S** | OMIQ e Cytobank (SVG, PDF) |
| Estilos de figura salvos e reaplicáveis | N? | **S** | — |
| Rampa de cor à escolha, inclusive cores próprias | P | **S** | — |
| Editor de layout com vários painéis | S | N | todos os comerciais |
| Relatório em lote | S | N | todos os comerciais |
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

8. **Testes entre grupos.** t, Mann-Whitney, ANOVA e Kruskal-Wallis sobre a frequência dos
   clusters e sobre as estatísticas da tabela, com correção para múltiplos testes e o valor de p
   no gráfico. O FlowJo não tem.
9. **Saída para o Prism.** Tabela já no formato de colunas agrupadas.
10. **Trabalhar os clusters.** Renomear, juntar, anotar; transformar um cluster ou um laço em
    população; sugerir nome pelo perfil de marcadores.
11. **Mais um método de cluster e um de redução.** PhenoGraph / Leiden, que não pede o número de
    clusters, e PCA.
12. **Normalização de lote.** Este experimento foi adquirido em três dias (9, 13 e 19 de novembro),
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
