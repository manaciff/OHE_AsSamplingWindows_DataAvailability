# Registro de correções, 24 de agosto de 2026

Sessão dedicada à revisão da lista de referências, à Figura S10 e à padronização
dos seis scripts. O manuscrito e o suplementar desta data foram aplicados sobre a
revisão de linguagem feita pela autora.

---

## 1. A revisão de linguagem da autora

**Verificação executada.** Extraí todos os tokens numéricos do corpo do
manuscrito, antes e depois da revisão, e comparei os dois conjuntos.

    643 tokens numéricos antes, 643 depois; nenhum valor removido, nenhum acrescentado.

Nenhum resultado foi alterado. Em seguida comparei sentença a sentença, sinalizando
apenas os pares em que mudou um número ou uma palavra de direção (excluded,
bounds, above, below, more, less, not, zero, positive, negative). Quatro
mudanças alteraram o sentido:

| Onde | Antes | Depois da revisão | Problema |
|---|---|---|---|
| 3.4 | "excluded zero in five of the six specifications" | "was excluded from five of the six specifications" | Inverte o sentido. "Excluded zero" diz que o intervalo de confiança não contém zero. "Was excluded from" diz que o preditor esteve ausente de cinco modelos, o que contradiz a Tabela S28a, onde a proximidade florestal está nos seis. |
| 3.4 | "a quantity the response bounds from above" | "a quantity that bounds the response from above" | Inverte a geometria. A floresta dentro de um envelope não pode exceder o envelope, de modo que é a resposta que limita a área de classe, e não o contrário. A Tabela S28e enuncia na direção certa, então os dois documentos passaram a discordar. |
| 4.2 | "although it carries no variation there that the other forest metrics do not" | "However, it shows no variation where the other forest metrics do not" | A afirmação é que a proximidade não carrega variação **única**, não que não carregue nenhuma. |
| 4.5 | "the most divided forest of the assemblage" | "the most fragmented forest of the assemblage" | O artigo reserva "fragmentation" para fragmentação per se, e a Seção 4.1 declara que este desenho não a mede. O uso solto convida à leitura que o artigo descarta. |

As quatro foram revertidas. Além disso, a revisão fundiu runs do documento e
desfez o itálico de cinco nomes de táxon no manuscrito e de um no suplementar;
todos foram restaurados.

No suplementar, duas mudanças de substância:

- "Tables S32 and S33 report the **recompilation** of the class-level metrics"
  voltou a "recomputation". As métricas foram recalculadas a partir dos rasters
  por unidade; recompilar seria reuni-las de novo, o que é outra afirmação.
- "The intercept is not constant across models because the predictors are
  standardized, but the fitted models differ" passou a "The intercept is not
  constant across models: the predictors are standardized, but the fitted models
  differ in composition". A pontuação anterior fazia da padronização a causa.

---

## 2. Referências

### 2.1 O protocolo do envelope é de Macedo et al. (2019)

A Seção 2.4 de Macedo et al. (2019) descreve o mesmo procedimento executado
aqui: mapa de densidade por kernel com `adehabitatHR`, protocolo ad hoc de Kie
(2013) partindo do h de referência e decrescendo 0,01, função `predict`,
*Panthera onca* na Mata Atlântica com 37.825 km² (Paviolo et al., 2016) como
referência para recuperar o h, o mesmo h aplicado a todos os táxons, recorte das
projeções fora da Mata Atlântica, de cidades e de corpos d'água, e projeção
Albers equivalente.

Até esta data o manuscrito citava Macedo et al. (2019) apenas como fonte dos
registros de ocorrência. A Seção 2.3 passa a declarar a filiação do protocolo, e
a Seção 2.3.2 marca a única diferença: o envelope é delimitado no isópleta de
99%, enquanto Macedo et al. retiveram a distribuição de utilização inteira.

Isso responde à pergunta sobre o bloco C da planilha. As referências de AOH e
AOO estão lá como **leitura de defesa**, não como candidatas a citação. O OHE não
deriva delas: deriva de Macedo et al. (2019), que chamam o resultado de "área
atualmente ocupada" e não usam nenhum dos dois termos. O parágrafo da Seção 2.3
existe para desambiguar um termo que a fonte nunca empregou, e um parágrafo de
desambiguação não precisa de seis referências.

### 2.2 Três citações acrescentadas

| Referência | Onde | Por quê |
|---|---|---|
| Keith et al. (2018) | 2.3 | Sustenta "No IUCN threshold is applied to any value reported here". Medidas de tamanho de área são dependentes de escala, e um limiar calibrado numa resolução não transfere para outra. Sem ela a frase era uma afirmação nua. |
| Püttker et al. (2020) | 4.1 | Fecha o parágrafo sobre a hipótese da quantidade de habitat. Modelos de equações estruturais sobre 1.097 sítios da Mata Atlântica encontram efeito negativo da fragmentação em níveis intermediários de quantidade de habitat, mediado por borda. O que falta a este desenho é a medida independente de quantidade, não a relevância da pergunta. |
| Magioli et al. (2015) | 4.5 | Contrapeso à segunda recomendação. A diversidade funcional de mamíferos da Mata Atlântica cresce com o tamanho da mancha até cerca de 2.050 ha, e a área média de mancha dentro dos envelopes descritos aqui vai de 12,6 a 132,2 ha, uma ordem de grandeza abaixo. Manchas pequenas bem posicionadas não substituem manchas grandes. |

As três sustentam um número que já estava no texto. As duas que aguardavam
decisão deixaram de estar destacadas em amarelo.

### 2.3 O que fazer com o resto do bloco C

Cinco entradas continuam como leitura, e uma está arquivada no lugar errado:

- **Marsh et al. (2023)**, efeito do esforço amostral sobre estimativas de área
  de distribuição, não é sobre AOH. É sobre viés amostral, que é o assunto do
  Apêndice S1. Deveria migrar para o bloco de viés amostral da planilha.
- **Palacio et al. (2021)** está descrita como "o protocolo mais próximo do que
  você executou". Não é: é um fluxo dedutivo de mapa de distribuição mascarado
  por habitat, a família oposta ao kernel indutivo sobre registros. Rebaixar.
- **Anderson et al. (2023)**, **Cazalis et al. (2024)** e **Alvarenga et al.
  (2025)** são contexto de arguição, não de texto.

---

## 3. Figura S10

`Figure_PIII_02_Coefficient_ForestPlot.tiff` era escrita pelo script 03 e não
aparecia em nenhum dos dois documentos. Ela mostra os seis coeficientes do
conjunto candidato a priori, com asteriscos de significância. Três dessas seis
barras, cobertura florestal, densidade de manchas e a interação, são exatamente
os coeficientes que o artigo conclui que não devem ser interpretados.

A figura entrou no suplementar como Figura S10, refeita a partir de
`Table_S_PartIII_FinalCoefficients.csv`, com o painel dividido em dois: os três
termos livres da resposta em azul cheio e com a marca de significância, os três
normalizados pela resposta em cinza vazado e sem marca. Uma estrela ao lado de um
coeficiente que o nulo de permutação reproduz é uma contradição impressa na
página. O script 03 foi alinhado para reproduzir essa versão.

---

## 4. Padronização dos scripts

### 4.1 Defeitos corrigidos

| Script | Defeito | Efeito |
|---|---|---|
| 00 | `stop(sprintf(...))` sem `paste0()`, de modo que a segunda string virava o primeiro argumento de formato | O guard de contagem de rasters disparava com `invalid format '%d'` em vez da sua mensagem |
| 00 | `canonical_taxon()` devolvia `NA` em silêncio para um nome desconhecido, enquanto as outras três cópias interrompem | Um nome fora da lista propagava para o join em vez de parar a execução |
| 02 | `vif_final` podia descrever um modelo diferente do conjunto retido quando o laço terminava por esgotar iterações ou por sobrar menos de dois preditores | O vetor obsoleto ia direto para `TableII_1_VIF_*.csv`. Na execução atual o laço saiu pelo ramo correto, então as tabelas exportadas estão certas |
| 02 | Comentário da Seção 9 afirmava que as frações de `varpart()` são lidas por posição | O código faz o contrário, deliberadamente, desde a correção de 22 de agosto. O leitor encontrava primeiro o comentário errado |
| 03 | Manifesto de exportação anunciava `Table_S_PartIII_OLS_Diagnostics.csv`, que nenhuma linha escreve, e descrevia o efeito aleatório como `(1 \| UA_ID)` quando o código usa `(1 \| GENUS)` | Refeito a partir do próprio código: 20 tabelas e 8 figuras, conferidas contra as chamadas de escrita |
| 05 | `MuMIn::AICc()` usado em seis lugares sem `MuMIn` no bloco de pacotes | Dependência não declarada, invisível para quem monta o ambiente pelo cabeçalho |
| 05 | Cabeçalho e `run_all.R` diziam "duas permutações de 9.999"; o segundo teste roda 2.000 | Documentação corrigida; o código não mudou |
| 99 | Cabeçalho dizia "não execute", mas o arquivo era executável e sobrescrevia três CSV de entrada da cadeia | Guarda acrescentada: interrompe a menos que `RUN_LEGACY_99 <- TRUE` seja definido antes |

### 4.2 Código morto removido

`01`: paletas `colors_locomotion` e `colors_diet` e a função `format_mean_sd()`,
nunca usadas; `library(ggtext)` e `library(patchwork)`, anexados sem uso.
`02`: `geo_path` e `shp_path`, mortos desde a remoção da reconstrução de
centroides; o bloco comentado da PCA e o bloco `if (FALSE)` da seleção
progressiva, este último referindo `data_z`, que a Seção 7 só cria abaixo dele, de
modo que ativá-lo falharia.
`03`: o operador `%||%`, definido e nunca usado; todo o ramo de erro-padrão
robusto HC3, inalcançável porque `use_robust_se` é sempre `FALSE`, e com ele
`library(sandwich)`; `interaction_term`, `use_glmm` e `use_robust_se`;
`library(stringr)`, que já vem no tidyverse.
`01`, `02`, `03`: vírgulas pendentes em `labs()` onde um título foi retirado, e
oito linhas de `theme()` que formatam `plot.title` e `plot.subtitle` em gráficos
que não têm nem título nem subtítulo.

### 4.3 Padronização

- `sessionInfo()` nu em `00`, `01` e `05` não imprimia nada sob
  `source(echo = FALSE)`, que é como `run_all.R` chama os scripts. Os seis passam
  a imprimir e a gravar `session_info_partNN.txt`.
- `04` era bilíngue. Os identificadores locais em português passaram a inglês.
  **Verificação:** executei o script renomeado sobre `Data_Raw_WithCoords.csv` e
  comparei a saída com a reprodução verificada de 23 de agosto. Todos os valores
  numéricos são idênticos; as únicas diferenças são os rótulos AOH que passaram a
  OHE de propósito.
- `05`: `THIN_RADIUS_M`, `N_THIN_REPS` e o teto do nulo de AICc, declarados no
  meio do arquivo, foram reunidos no bloco de constantes do topo.
- `03`: `N_UNITS_EXPECTED` foi do meio do arquivo para o bloco de constantes;
  quatro consultas cruas a `predictor_labels[...]` passaram a usar
  `label_predictors()`, que é o helper escrito para evitar `NA` silencioso.
- `run_all.R`: a mensagem final afirmava "Every analysis used 67 sampling units".
  Passa a ler a contagem do arquivo que a cadeia acabou de escrever.

### 4.4 Seção 16 do script 03, removida

O arquivo carregava 117 linhas de notas interpretativas escritas para uma versão
anterior do estudo. Elas contradiziam o próprio script: chamavam a resposta de
AOH; afirmavam que as 67 observações se distribuem por **12** unidades amostrais
e que o efeito aleatório é a unidade amostral, quando são 67 unidades e o efeito
aleatório é o gênero; citavam 52,9% de variância não explicada, "~47%" e "~31%",
VIF máximo 2,2 e acurácia do MapBiomas de 87%, todos superados; e remetiam a
seções 6.5, 8.5, 8.6 e 10, que não existem no arquivo. Notas que contradizem o
código em que estão são piores do que nenhuma nota.

Foram substituídas por um bloco curto que aponta onde cada decisão está
registrada hoje: Seções 2.5.1, 2.5.5 e 4.1 do manuscrito, `docs/TABLE_MAP.md`,
`docs/SCRIPTS.md` e os registros de correção. O texto retirado está arquivado
fora do repositório, com este registro.

### 4.5 O que deixei como está, e por quê

- **`canonical_taxon()` em quatro cópias.** Consolidar num arquivo comum muda a
  estrutura da cadeia. Alinhei o comportamento das quatro; a duplicação continua
  e está anotada.
- **Refits repetidos no script 05.** O mesmo modelo é ajustado cinco vezes.
  Custa segundos e mexer em código de ajuste arrisca os resultados.
- **Identificadores em português no script 03.** São treze, espalhados por 2.200
  linhas. O ganho é de estilo e o risco não é.

---

## 5. Pendências

1. **Rodar o script 04.** `TableS28_Geometric_Coupling.csv` na sua pasta ainda é
   de 23/08 às 12h07, produzido pela versão antiga, e traz 0,0290 e o rótulo
   `log(AOH)`. A Tabela S26 do suplementar já traz 0,0277. São segundos de
   execução, mas até lá o arquivo do repositório discorda do documento.
2. Regerar `Figure_PIII_02_Coefficient_ForestPlot.tiff` na próxima execução do
   script 03, para que o arquivo do repositório seja a versão dividida em dois
   painéis que está no suplementar como Figura S10.
3. Preencher os três marcadores do `CITATION.cff` e as duas afiliações.
