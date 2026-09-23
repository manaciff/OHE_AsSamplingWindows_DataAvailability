# Registro de correções da apostila técnica

**Objeto:** `Apostila_Tecnica_Analises_TCC_v3.docx`, de 24 de agosto de 2026
**Resultado:** `Apostila_Tecnica_Analises_TCC_v4.docx`, de 17 de setembro de 2026
**Confrontada com:** os 69 arquivos de `Outputs/Manuscrito/`, o `Manuscript_Revisado_15.09.docx`, o `Supplementary_Material_REVISED_15.09.docx` e os scripts em `Repositorio_GitHub/R/`

A versão 3 citava números das saídas de 30 de julho. Desde então a cadeia foi rodada de novo e o manuscrito foi revisto duas vezes. Trinta e sete divergências foram encontradas. Onze delas são materiais: uma banca que abrisse o suplementar encontraria contradição direta com o que a apostila afirmava.

---

## A. Divergências materiais

### A1. Os três números-âncora da Seção 10 estavam todos errados

A apostila mandava conferir, ao fim da reprodução, o R² ajustado da RDA, a fração de paisagem da partição e o ICC do gênero. Os três valores impressos eram de uma rodada anterior.

| Âncora | Apostila v3 | Saída atual | Arquivo |
|---|---|---|---|
| R² ajustado da RDA | 0,584 | **0,905** | `PartII/TableII_4_RDA_GlobalSummary.csv` |
| Fração de paisagem pura | 62,6% | **51,6%** | `PartII/TableII_3_VP_Fractions.csv` |
| ICC do gênero | 0,138 | **0,0758** | `Table_S_PartIII_DependenceDiagnostic.csv` |

Quem seguisse o roteiro concluiria que a cadeia não reproduz. Os três valores corretos já apareciam nas Seções 6.3 e 6.4 da própria apostila, o que tornava a contradição visível dentro do documento.

### A2. A Seção 10 começava no item 6 e citava scripts que não existem

O capítulo abria em "6. Confirmar que Data_Raw_FINAL.csv...", sem os itens 1 a 5, e mandava rodar `Script_Part_I_v9.R`, `Script_Part_II_v4_CORRIGIDO.R` e `Script_Part_III_v2_CORRIGIDO_v3.R`. Nenhum dos três existe. A cadeia atual tem seis scripts numerados, de `00` a `05`, chamados por `run_all.R`, e é o que a Seção 5 da própria apostila descreve. O capítulo foi reescrito inteiro, em seis passos, com os arquivos e os números conferidos.

### A3. O betadisper foi lido ao contrário

**Seção 4.9 da v3:** "a dispersão foi homogênea para gênero, espécie e locomoção e heterogênea para dieta (p = 0,004)".

**`PartI/TableS3_Betadisper.csv`:** heterogênea em gênero (p = 0,0424), em espécie (p = 0,0405) e em dieta (p = 0,0044); homogênea apenas em locomoção (p = 0,1261).

A afirmação da Seção 4.9 contradizia a Seção 6.0 do mesmo documento, que já trazia a leitura correta. A consequência é de interpretação: a PERMANOVA de gênero e de espécie responde à posição dos centroides e também à dispersão interna dos grupos, e é assim que deve ser apresentada. Corrigido em 4.9 e ressalvado em 6.2.

### A4. A amplitude da resposta era a da coluna arredondada

**Seção 4.5 da v3:** "varia de 1.910 a 387.619 ha", e "a média é 10,38 e o desvio é 0,95".

**`Table_S_PartIII_DescriptiveStatistics.csv`:** 1.906,73 a 386.509,44 ha; log(OHE) com média 10,3735 e desvio 0,9392.

A origem do erro é informativa e entrou na apostila: `Data_Raw_FINAL.csv` guarda `UA_Area_ha` arredondada ao hectare, de 1.910 a 387.619 ha, com média de log de 10,3784 e desvio de 0,9383. O script 01 substitui essa coluna pela área de precisão total de `AOH_area_full_precision.csv`, e é essa que entra nos modelos. A v3 citava a coluna errada em 4.5 e a certa na tabela de 6.1. A v4 explica a substituição na Seção 5 e usa a diferença entre as duas colunas como conferência no passo 4 da Seção 10.

### A5. O erro de validação cruzada estava desatualizado em dois lugares

| Onde | v3 | Saída atual | Arquivo |
|---|---|---|---|
| 4.13, Gama contra MQO | 35.189 e 42.040 ha | **39.977 e 48.572 ha** | `Table_S_PartIII_ModelComparison.csv` |
| 4.21, acoplado contra reportado | 40.102 e 128.197 ha | **39.977 e 127.921 ha** | `PartV/TableS_PartV_CV_And_BIC.csv` |

A v3 também dizia que o Gama era "20% melhor". A redução é de 17,7% no leave-one-out e de 18,5% no leave-one-genus-out; dito ao contrário, o erro do MQO é 21,5% maior. As duas formulações são verdadeiras e dizem coisas diferentes. A v4 usa "cerca de 18%, ou, dito ao contrário, cerca de 21% maior", e avisa para não misturar as duas.

### A6. Percentil e valor de p foram trocados no teste de permutação

**Seção 4.16 da v3:** os coeficientes "caíram todos dentro da distribuição nula, nos percentis 2,9, 7,5 e 15,5".

**`PartV/TableS_PartV_NullSimulation.csv`:** os percentis são 2,8, 6,9 e 10,4. Os números que a v3 citava, 2,9, 7,5 e 15,5, são os valores empíricos de p bilaterais, 0,0294, 0,0748 e 0,1545, expressos em porcentagem.

Não é apenas um erro de valor: é a troca de duas quantidades distintas. O manuscrito de 15.09 usa os percentis 2,8, 6,9 e 10,4, de modo que a apostila dizia uma coisa e o artigo, outra. A v4 traz os dois conjuntos, nomeados.

### A7. Os coeficientes do modelo a priori estavam arredondados para baixo

**Seção 4.16 da v3:** PLAND −1,302, PD −0,642, interação −0,712.
**`Table_S_PartIII_ConfirmatoryModel.csv`:** −1,3062, −0,6440, −0,7138, ou seja, **−1,306, −0,644 e −0,714**.

A Tabela 3 do manuscrito traz os valores corretos. Também nesta seção: o erro máximo da identidade algébrica é de **0,028**, e não 0,029 (`TableS28_Geometric_Coupling.csv`); a vantagem do modelo acoplado é de **17,89** unidades de AICc, alcançada em **38,1%** das permutações, e não 17,8 e 38,8% (`PartV/TableS_PartV_Decision.csv`, Q5); e a correlação bruta de PLAND com log(OHE) é de **0,129**, e não 0,130.

### A8. A tabela de pesos do índice de conservação omitia um dos nove táxons

A tabela da Seção 6.8 listava oito táxons, cujos pesos somam 18, e o texto logo abaixo falava num "máximo teórico de 21". A inconsistência tinha causa: **Tayassu pecari**, peso 3, não estava na tabela.

A Seção 5.2 do suplementar atribui peso 3 a *Brachyteles arachnoides*, *Brachyteles hypoxanthus*, *Bradypus crinitus* e *Tayassu pecari*; peso 2 a *Alouatta guariba*, *Leopardus wiedii*, *Mazama* e *Puma concolor*; e peso 1 a *Myrmecophaga tridactyla*. Soma 21. A linha foi devolvida, e o nome da preguiça foi completado para *Bradypus crinitus*, como no suplementar.

### A9. O β da matriz agropecuária misturava dois modelos

**Seção 6.6 da v3:** "com β de 0,51 ... Todos os VIFs ficaram abaixo de 2,3", seguido de uma tabela em que o mesmo preditor aparece com β de 0,400 e VIF de 1,96.

São dois ajustes diferentes. O β de 0,508, com VIF de até 19,3, é do melhor subconjunto de cinco preditores do conjunto a priori, que contém métricas normalizadas pela resposta (`Table_S_PartIII_Sensitivity_MatrixModel.csv`, Tabela S23). O β de 0,400, com VIF de 1,96, é do modelo livre de acoplamento (`PartV/TableS_PartV_Matrix_Uncoupled.csv`, Tabela S40), que é o que o manuscrito reporta e o que a tabela da apostila exibia. Citar os dois lado a lado é um erro verificável. A v4 usa o 0,400, informa o 0,508 como pertencente a outro modelo, e diz qual é qual.

### A10. A Seção 4.4 afirmava que todos os VIFs ficam abaixo de 5

A tabela de VIF da v3 trazia quatro linhas, e as quatro estavam erradas. Pior: a conclusão "todos abaixo de 5" contradiz o resultado metodológico central do trabalho, exposto em 4.16, de que a cobertura florestal entra no modelo a priori com VIF de 16,8.

| Conjunto | v3 | Saída atual | Arquivo |
|---|---|---|---|
| RDA | 3,88, dimensão fractal | **3,20**, dimensão fractal | `PartII/TableII_1_VIF_RDA.csv` |
| Partição de variação | 3,56, distância ao vizinho | **3,92**, dimensão fractal | `PartII/TableII_1_VIF_VP.csv` |
| Modelo univariado | 2,17 | **16,82**, cobertura florestal | `Table_S_PartIII_VIF.csv` |
| Modelo de sensibilidade | 2,26 | **1,96** | `PartV/TableS_PartV_Matrix_Uncoupled.csv` |

A v4 traz os quatro valores e usa o contraste entre 16,82 e 1,96 como argumento, que é o que ele é.

### A11. A Seção 7 remetia a uma numeração de tabelas que não existe mais

O guia listava S1 a S27 e atribuía a S26 e S27 os pesos e a distribuição do índice de conservação. No suplementar atual, o índice está em **S36 e S37**, a **S27** é a contagem de unidades por análise e a **S16** compara as seis especificações, não mínimos quadrados contra modelo misto. Faltavam ainda vinte tabelas, entre elas as S26 e S28a a S28f, que documentam o acoplamento geométrico, isto é, justamente o achado central. A seção foi refeita, com as 48 tabelas, de S1 a S42b, e as dez figuras.

A legenda da Figura S5 também estava errada em três valores: o R² anotado no painel é **0,662** na escala de hectares e **0,557** na logarítmica, não 0,13, 0,43 e 0,54, e a maior unidade tem **386.509 ha** (`Table_S_PartIII_ObsPred_R2.csv` e legenda do suplementar).

---

## B. Contradições internas

| # | Onde | O que dizia | Correção |
|---|---|---|---|
| B1 | 1.3 e 2.2 contra 3.0 e 6.1 | "67 fragmentos, dos quais 65 foram incluídos"; "Data_Raw_FINAL.csv, 65 linhas" | 67 unidades em todas as análises; o arquivo tem 67 linhas, conferido |
| B2 | 4.14 contra 6.5 | "só um modelo entrou no conjunto competitivo" | Dois modelos, separados por 1,51 unidade de AICc (`Table_S_PartIII_CompetitiveModels.csv`) |
| B3 | 4.14 e 4.19 | "16 subconjuntos dos quatro preditores" | 32 subconjuntos de cinco preditores mais a interação (`Table_S_PartIII_Dredge_FullTable.csv`, 32 linhas) |
| B4 | 4.13 contra 6.4 | ICC de 7,5% num parágrafo, 0,0758 no outro, e uma caixa explicando o ICC com o valor antigo de 0,138 | 7,6% nos três lugares |
| B5 | 6.7 | "o efeito mais forte em um gênero isolado é o de *Brachyteles*", e três linhas abaixo "o β de 2,29 em *Myrmecophaga*, o maior da tabela" | *Brachyteles* é o único apoiado; o β de *Myrmecophaga* é 2,26, o maior e não significativo |
| B6 | 9, glossário contra 4.18 | "AOH, Área de Habitat" e "unidade amostral: há 65 delas" | "OHE, envelope de habitat ocupado", com a distinção em relação a Brooks et al. (2019) e à AOO da IUCN; 67 unidades |

---

## C. Valores pontuais corrigidos

| Seção | v3 | Correto | Fonte |
|---|---|---|---|
| 4.8 | estresse do NMDS 0,089 | **0,084** | `PartI/TableS_NMDS_Summary.csv` |
| 4.8, caixa | densidade de borda, R² de 0,144 | **0,165** | `PartI/TableS2_Envfit_Configuration.csv` |
| 4.13 | jacobiano 695,35 | **695,03** | soma de log(OHE) em `Data_Raw_WithCoords.csv`; 695,35 é a soma da coluna arredondada |
| 4.17 | proximidade "0,741 vira 0,747" | **0,742 vira 0,748** | `PartV/TableS_PartV_Proximity_Conditioned.csv` |
| 4.20 | deslocamento máximo "mais de 0,16"; isolamento p = 0,223; forma p = 0,017 | **0,16 exatos**; **0,224**; **0,016** | `Table_S_PartIII_SpatialSensitivity.csv` |
| 4.21 | Tabela S30c | **Tabela S28c** | suplementar |
| 4.23 | índice sustenta a Figura 7 | **Figura 8** | suplementar, Seção 5.4 |
| 6.2 | agropecuária vai "a 68,6% em *Alouatta*" | *Mazama* tem a mediana mais alta, **70,1%** | `PartI/Table3_Medians_Composition_PLAND.csv` |
| 6.6 | proximidade "de 1,43 para 1,48 sem a forma" | parte de **1,36** nos dois casos | `Table_S_PartIII_StabilityCheck.csv` |
| 6.6 | Spearman PD × FRAC de 0,723 | **0,729** | idem |
| 6.7 | *Mazama*: PD −0,57 e ENN −1,17 | **−0,58** e **−1,15** | `Table_S_PartIII_EffectByGenus.csv` |
| 6.5 | tabela de importância com quatro preditores | seis, incluindo PLAND (0,9996) e a interação (0,9987) | `Table_S_PartIII_VariableImportance.csv` |
| 8 | "táxons com massa superior a 1,5 kg" | **1,0 kg**, como a Seção 2.2 do manuscrito | `Manuscript_Revisado_15.09.docx` |

---

## D. Estrutura

1. **A Seção 4.15 estava vazia.** O título "Três verificações que não estavam obrigadas a existir" não tinha corpo, e o texto correspondente, com a tabela das três verificações, estava solto ao fim da Seção 4.18, que trata do nome da resposta. O conteúdo voltou para o seu lugar.
2. **Parágrafos duplicados.** A Seção 3.2 trazia duas vezes "a correção não usou a tabela auxiliar"; a Seção 3.3, duas vezes "a consequência é metodológica"; a Seção 3.0, o mesmo argumento do quadro que fecha a Seção 3.3. As repetições foram fundidas.
3. **A Seção 4.13 estava fora de ordem.** Abria pelo resultado, voltava à alternativa descartada e terminava com "o que foi feito" descrevendo o procedimento antigo, que a própria seção diz ter sido abandonado. Refeita na ordem do resto do capítulo. Foi acrescentado um alerta que faltava: o modelo retido na comparação das seis especificações não é o modelo reportado no artigo, e os 55,1% de desvio explicado pertencem a ele, não à Tabela 4.
4. **A Seção 4.12 abria por uma observação de rodapé**, antes do "o que foi feito". A observação foi para o fim.
5. **A Seção 8 tinha a frase de fecho no meio da lista.** Foi para o fim.
6. **A tabela de scripts da Seção 5** dava `Data_Raw_FINAL.csv` como entrada da Parte III. O script 03 recusa esse arquivo explicitamente e exige `Data_Raw_WithCoords.csv`. A tabela foi refeita, com as seis partes.
7. **Uma frase invertia o sentido.** Em 4.13: "Isso tem dois defeitos. Nunca compare as três famílias em pé de igualdade e deixe que um diagnóstico marginal decida uma questão estrutural." O imperativo mandava fazer o que a frase queria condenar. Corrigido para "nunca compara... e deixa...".
8. **Pontuação em 6.8:** "Segundo o índice, mede-se a coocorrência" e "Quarto é uma camada de triagem", nos dois casos com a vírgula no lugar errado.

---

## E. Conteúdo removido

- **Seção 4.23**, sobre o rótulo NA na Figura PIII_03, cerca de 330 palavras. É depuração de código, não material de arguição. A Seção 4.24 passou a ser 4.23.
- **Seção 5.2**, com a tabela das quatro correções de código para modelos `lmer`, e o quadro sobre predição populacional, cerca de 380 palavras. O quadro descrevia uma escolha relativa a modelo misto; o modelo da Figura S5 é o Gama acoplado, sem efeito aleatório, de modo que o quadro já não se aplicava.
- **Item obsoleto no quadro de pendências:** "o critério de exclusão das duas unidades amostrais". Nenhuma unidade é excluída: as duas foram recuperadas pelo script 00, conforme a Tabela S33.
- Trechos repetidos entre 4.11 e 6.3, entre 4.12 e 6.3, e entre 4.13 e 6.4 foram reduzidos a remissões.

---

## F. O que permanece sem verificação

Três valores citados não têm contrapartida em arquivo exportado. Nenhum altera as conclusões, e os três estão sinalizados na própria apostila.

1. **VIF de 982 da cobertura florestal** na primeira iteração da eliminação por VIF, Seção 4.12. É impresso no console do script 02.
2. **Escore de biplot de 0,68** da proximidade no primeiro eixo da RDA, Seção 6.3.
3. **Correlação de 0,94** entre o envelope e a área de floresta, Seções 4.17 e 4.18. O valor consta da Tabela S30 do suplementar, na coluna de relação com a resposta, mas não como número exportado pela cadeia em R.

---

## G. Dois pontos a confirmar antes da defesa

1. **O limiar de porte.** O manuscrito de 15.09 define porte médio a partir de 1,0 kg, citando Chiarello (2000); até 2 de setembro o limiar era 1,5 kg. A apostila foi alinhada ao manuscrito. Confirme que Chiarello (2000) sustenta 1,0 kg, porque é a única fonte dada.
2. **O desvio explicado do modelo reportado.** A apostila diz "cerca de 36%". O valor 0,3609 está exportado para a especificação "S4 mais densidade de manchas da matriz", cujo preditor adicional é nulo (β = −0,024; p = 0,808), de modo que a deviância é praticamente a do S4. O desvio do S4 isolado não está exportado em nenhum arquivo. Se a banca pedir o número exato, ele precisa ser lido no script.
