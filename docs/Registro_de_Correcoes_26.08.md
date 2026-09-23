# Registro de correções, 26 de agosto de 2026

Trabalho feito em resposta ao parecer recebido e a uma revisão independente da
metodologia. Arquivos entregues:

- `Manuscript_REVISED_26.08.docx`
- `Supplementary_Material_REVISED_26.08.docx`
- `R/08_sampling_effort_sensitivity.R`
- `Outputs/Manuscrito/PartVI/` com a contagem por unidade e as Tabelas S42 e S42b
- `Resposta_ao_Parecer_26-08-2026.md`

Cada edição abaixo traz o texto original, o texto revisado e o motivo. Nenhum
valor de nenhuma tabela existente foi alterado, porque nenhum estava errado.

---

## 1. O que o parecer acertou e o que não acertou

O parecer levanta três diagnósticos estatísticos e duas questões de escopo. Os
três diagnósticos descrevem problemas que o manuscrito já identifica, mede e
resolve, citando a mesma literatura que o parecer invoca. A crítica se apoia na
Tabela S18, que é tabela de diagnóstico do modelo que o manuscrito
declaradamente não interpreta. A legenda dessa tabela, porém, chamava aquele
modelo de "the retained model", rótulo herdado de versão anterior, e foi isso
que induziu a leitura. O defeito é do documento, e foi corrigido.

Das cinco críticas, uma procede inteiramente: as limitações do índice de valor de
conservação estavam apenas no suplemento, enquanto o texto principal usava o
índice sem nenhuma delas. Foi corrigido.

A recomendação de mudar a narrativa foi recusada. O objetivo declarado é
descrever composição e configuração, a ênfase metodológica que o parecer pede já
existe na Seção 4.1 e no primeiro parágrafo da Conclusão, e pedir outra
apresentação é pedir outro artigo.

---

## 2. Achado próprio: dependência entre número de registros e área do kernel

Ao rever a metodologia identificou-se uma dependência que o parecer não menciona
e que é mais séria do que as três apontadas. A resposta é a área de uma
distribuição de uso estimada por kernel na isópleta de 99% com parâmetro de
suavização único. A parâmetro fixo, essa área cresce com o número e a dispersão
dos registros que a produziram, e o esforço em registros secundários concentra-se
perto de estradas, estações de pesquisa e áreas protegidas, que são também onde a
floresta persiste. O Apêndice S1 tratava da irregularidade espacial do esforço,
que é outro mecanismo.

Contaram-se os registros retidos que caem dentro de cada unidade, 503 dos 606, e
o modelo da Tabela 3 foi refeito com o logaritmo dessa contagem como covariável.

| Preditor | β, Tabela 3 | β, com esforço | IC 95% | Razão por DP |
|---|---|---|---|---|
| Proximidade | 0,742 | 0,412 | de 0,205 a 0,619 | 1,51 |
| Complexidade de forma | 0,602 | 0,489 | de 0,313 a 0,665 | 1,63 |
| Isolamento | 0,122 | 0,106 | de -0,045 a 0,256 | 1,11 |
| Esforço, log(registros + 1) | --- | 0,488 | de 0,331 a 0,646 | 1,63 |

As duas associações reportadas sobrevivem, com sinal preservado e intervalo
excluindo zero. A deviance explicada sobe de 36,1% para 57,0% e o fator de
inflação máximo de 1,93 para 2,32. O jackknife sobre as 67 unidades mantém a
proximidade entre 0,364 e 0,448 e a forma entre 0,460 e 0,523.

Dois fatos ficam declarados no texto. O primeiro é que log(registros + 1)
correlaciona 0,588 com log(OHE) e é o correlato marginal mais forte da resposta,
acima da proximidade, com 0,350. O segundo é que o coeficiente de proximidade cai
44% sob o controle.

A leitura adotada é que os dois conjuntos de coeficientes limitam as associações
em vez de um corrigir o outro. A contagem de registros é em parte mediadora da
mesma quantidade latente que a resposta mede, porque a premissa da Seção 2.3 é
que populações maiores mantêm áreas de uso mais amplas e populações assim geram
mais registros. Condicionar nela remove esforço e sinal ecológico ao mesmo tempo.
Os valores condicionais são piso, os da Tabela 3 são teto, e ambos excluem zero.

Apenas o subconjunto do Rio de Janeiro da compilação está arquivado na pasta,
enquanto os kernels vieram da compilação de bioma inteiro de Macedo et al.
(2019). A contagem mede densidade de registros dentro de cada envelope, não o
insumo completo do kernel. Dez unidades não contêm registro desse subconjunto,
razão pela qual a covariável é log(registros + 1) e a análise se repete nas 57
unidades com ao menos um.

---

## 3. Edições no manuscrito

### 1. Section 3.4, last sentence of the Table 3 paragraph

> **Antes:** the absolute class area of forest, a quantity that bounds the response from above
>
> **Depois:** the absolute class area of forest, a quantity the response bounds from above
>
> Motivo: The forest inside an envelope cannot exceed the envelope, so the response bounds the class area and not the reverse. Tables S28a, S28e and S28f state the correct direction, so the two documents currently disagree.

### 2. Section 3.4, new paragraph after the Table 3 paragraph

> **Antes:** (new paragraph)
>
> **Depois:** Both associations also survive a control for sampling effort. The number of retained occurrence records inside an envelope is the strongest single correlate of its extent in these data, at r = 0.59 on the logarithmic scale, against 0.35 for forest proximity and 0.13 for forest cover. Adding the logarithm of that count to the model leaves both coefficients positive with intervals that exclude zero, at 0.412 (0.205 to 0.619) for clustering and 0.489 (0.313 to 0.665) for shape complexity, so the envelope is 1.51 and 1.63 times larger per standard deviation rather than 2.10 and 1.83 (Table S42). The record count is in part a consequence of the same population extent that the response measures, so the two sets of coefficients bound the associations from above and from below rather than one correcting the other.
>
> Motivo: The response is a kernel area at fixed bandwidth, so it grows with the number of records that produced it. The new analysis in Table S42 shows both associations survive that control, at about two thirds of their unconditional magnitude.

### 3. Section 2.5.5, after the sentence on selection uncertainty

> **Antes:** (sentence added)
>
> **Depois:** The specification reported in Table 3 is the member of that set from which every metric normalized by the area of the sampling window has been removed. That criterion follows from the formulas of the metrics, and it can be applied before any model is fitted, so those estimates are not conditioned on a selection step either. The permutation test in Section 2.5.1 measures whether the restriction is necessary; it does not choose the predictors.
>
> Motivo: Table 3 is specification S4, declared a priori in the analysis script alongside S1 to S6 and defined by a rule that reads the metric formulas rather than the data. Saying so removes the reading that the permutation diagnostic selected the model, which is the post-selection objection.

### 4. Section 4.6, new paragraph after the first

> **Antes:** (new paragraph)
>
> **Depois:** A related dependence follows from the response itself. The envelope is a kernel utilization distribution taken at a fixed smoothing parameter, so its area grows with the number and the spread of the records that produced it, and effort in secondary records concentrates near roads, research stations, and protected areas, which are also where forest persists. Table S42 measures that dependence. Both reported associations keep their sign and exclude zero once the number of records inside the envelope is controlled for, at roughly two thirds of their unconditional magnitude, and because the record count is itself in part a consequence of population extent, the conditional estimates are a lower bound rather than a correction. Only the Rio de Janeiro subset of the compilation is archived here, so the count describes the density of records inside each envelope rather than the complete input to the kernel.
>
> Motivo: The dependence between record count and kernel area is not covered by Appendix S1, which addresses the spatial unevenness of effort. Stating it, with the quantity attached, closes the objection before a referee raises it.

### 5. Section 4.5, new paragraph before the numbered recommendations

> **Antes:** (new paragraph)
>
> **Depois:** One qualification applies to the Weighted Conservation Value Index used below. Its weights are ordinal and expert-derived, so the ranking is conditional on that judgment. It measures the weighted co-occurrence of habitat envelopes rather than the condition of the forest within them, so a degraded or defaunated remnant can rank highly on biogeographic position alone. It contains no measure of connectivity or forest cover, and it carries no information on land tenure, opportunity cost, restoration feasibility, or governance. The index is therefore a biological screening layer that shortlists candidate areas for further assessment, in combination with socioeconomic data and field validation, and not an allocation of conservation effort. Section 5 of the supplementary material sets out these limits in full.
>
> Motivo: Supplementary Section 5.5 states these four limits honestly, but the main text uses the index to justify the first recommendation and Figure 8 without any of them. This is the one point of the review that stands on its own.

### 6. Abstract

> **Antes:** (sentence added after the two headline ratios)
>
> **Depois:** Both persisted, at 1.51 and 1.63 times, when the number of occurrence records inside the envelope was controlled for.
>
> Motivo: The abstract announces 2.10 and 1.83. The conditional values belong beside them, so the headline is not read as unconditional.

### 7. Reference list, two paragraphs each holding two entries

> **Antes:** "Cushman SA, McGarigal K, Neel MC (2008). Parsimony in landsc ..." + "da Silva FR, Oliveira-Silva AE, Antonelli A, Carnaval AC, Pr ..."; "Legendre P, Gallagher ED (2001). Ecologically meaningful tra ..." + "Legendre P, Legendre LFJ (2012). Numerical ecology, 3rd edn. ..."
>
> **Depois:** each pair separated into two entries; author list corrected to "Oliveira-Silva AE, Antonelli A"
>
> Motivo: Two pairs of references sit in a single paragraph each, joined by a stray U+2029 character, so the reference list runs two entries short in any automatic check and the hanging indent does not apply to the second of each pair. The author list of da Silva et al. also carried a stray "de".

### 12. Section 3.4, deviance explained by the reported model

> **Antes:** the model finally reported in Table 3 ... accounts for 36.0%.
>
> **Depois:** ... accounts for 36.1%.
>
> Motivo: The reported model explains 0.360784 of the deviance, which prints as 36.1 percent, not 36.0. The figure for the coupled model in the same sentence, 55.1 percent, is correct at 0.550876, so the two were not rounded by the same rule. Table S42 now cites the same quantity, so the two documents would otherwise disagree.

---

## 4. Edições no material suplementar

### 8. Caption of Table S18

> **Antes:** Variance Inflation Factors for the terms of the retained model of log(OHE), computed with the car package (Fox & Weisberg, 2019).
>
> **Depois:** Variance inflation factors for the terms of the full a priori model of log(OHE), the model reported in Table S18a, computed with the car package (Fox & Weisberg, 2019). This is the coupled specification, which Table 3 of the manuscript supersedes for interpretation. The model reported in Table 3 contains no metric normalized by the response and has a maximum factor of 1.93, given as specification S4 in Table S28a.
>
> Motivo: The caption calls this the retained model. It is not: it is the coupled a priori model, whose estimates Table S18a reports and Table 3 supersedes. The wording invites a reader to attribute a factor of 16.8 to the model the manuscript actually reports, which has a factor of 1.93.

### 9. Legend beginning "Predictors were standardized before forming the product term"

> **Antes:** printed after Table S18b, the model-averaging table
>
> **Depois:** moved to immediately after Table S18, the variance inflation table
>
> Motivo: The legend explains the variance inflation factors, the threshold of 5 and the correlation of 0.66 between forest cover and forest proximity. Printed under Table S18b it reads as a legend for the model averaging, which reports no factor at all.

### 10. New Tables S42 and S42b, after the legend of Table S41

> **Antes:** (new)
>
> **Depois:** Table S42, sampling effort and envelope extent; Table S42b, robustness of that control
>
> Motivo: The dependence between the number of records and the area of a kernel at fixed bandwidth is not covered anywhere in the document. Appendix S1 treats the spatial unevenness of effort, which is a different matter. Both reported associations survive the control, so the analysis strengthens the paper and closes an objection that a referee would otherwise raise.

### 11. Appendix S1, end of the paragraph on sampling bias

> **Antes:** This limitation is revisited in Section 4.6.
>
> **Depois:** This limitation is revisited in Section 4.6. A second and distinct dependence arises from the response rather than from the predictors: at a fixed smoothing parameter the area of the kernel grows with the number and the spread of the records that produced it. Table S42 measures that dependence and reports the two associations of Table 3 of the manuscript with the number of records inside the envelope entered as a covariate.
>
> Motivo: Appendix S1 discusses the spatial unevenness of effort but not the dependence between the number of records and the area of the kernel, which is a different mechanism. The cross-reference sends the reader to the table that now measures it.

---

## 5. Edições nos scripts

### 12. `R/03_univariate_model_selection.R`, bloco de estatísticas descritivas

Os rótulos exportados em `Table_S_PartIII_DescriptiveStatistics.csv` eram
`AOH_ha` e `log_AOH`, enquanto a Tabela S15 do suplemento imprime `OHE_ha` e
`log_OHE`. Acrescentou-se um `mutate` que renomeia apenas os rótulos impressos.
Os nomes internos de coluna não mudam, portanto nenhum outro trecho da cadeia é
afetado. Pendência registrada em 25 de agosto.

### 13. `R/05_ratio_artefact_test.R`, duas ocorrências

A coluna `AOH_ratio_per_SD` passou a `OHE_ratio_per_SD`, em acordo com
`TableS_PartV_Matrix_Uncoupled.csv`, escrita pelo mesmo script. Duas ocorrências,
linhas 849 e 959. Nenhum outro arquivo lê essa coluna pelo nome, o que foi
verificado por busca em todo o repositório. Pendência registrada em 25 de agosto.

### 14. `R/08_sampling_effort_sensitivity.R`, novo

Produz as Tabelas S42 e S42b e a contagem de registros por unidade. Semente fixa
em 123, como o restante da cadeia. Lê `Coordenates_sp_Data.xlsx`,
`Data_Raw_WithCoords.csv` e os nove shapefiles de
`08_Dados_Especies/Dados_geo_especies/Sp_data_singlepart/`.

O casamento entre polígono e `UA_ID` é feito por área planar na projeção Albers
dos polígonos, que concorda com a área geodésica armazenada na tabela de análise
com erro relativo abaixo de 1,2% em todas as 67 unidades. O script exige que a
atribuição seja bijetiva dentro de cada táxon e que o erro relativo máximo fique
abaixo de 0,02, e para em vez de prosseguir sobre um casamento errado.

**Ainda não rodado no R.** Os números das Tabelas S42 e S42b vieram de uma
implementação independente em Python, que reproduziu a Tabela 3 até a terceira
casa decimal antes de qualquer coisa nova ser calculada: proximidade 0,7413
contra 0,7420, forma 0,6015 contra 0,6023, isolamento 0,1205 contra 0,1217, e
deviance explicada 0,3608 contra 0,3609. Rode o script uma vez no seu R para
confirmar antes de citar os valores em banca. Ele requer `sf` e `readxl`, e
exporta também um teste de resíduos por simulação com `DHARMa` que não consta da
Tabela S42b; se quiser incluí-lo, cole a linha do CSV gerado.

O script não foi acrescentado a `run_all.R`, pelo mesmo motivo que os scripts 06
e 07 não estão: depende de pacotes fora da cadeia inferencial. Acrescente-o
depois de confirmar que roda.

---

## 6. O que continua pendente

Herdado do registro de 25 de agosto e não resolvido aqui:

1. Depositar os 67 rasters por unidade, 2,6 GB, em repositório com DOI, e
   registrar o DOI no `CITATION.cff`.
2. Preencher `date-released`, `repository-code` e `doi` no `CITATION.cff`.
3. Decidir as afiliações dos coautores 2 e 3, hoje marcadas em amarelo.

Novo:

4. Rodar `R/08_sampling_effort_sensitivity.R` no R e conferir as Tabelas S42 e
   S42b contra os valores impressos no suplemento.
5. Decidir se a Figura 6 e a Figura S10 ganham o par de coeficientes
   condicionais ao esforço, hoje presentes apenas na Tabela S42.
