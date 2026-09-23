# Registro das correções aplicadas em 22 de agosto de 2026

Documento de acompanhamento do parecer de auditoria de 22/08/2026. Lista o que foi alterado, onde, e com que evidência. Ao final, o que ainda depende de você.

Todos os valores novos foram obtidos por execução independente em R 4.4 com *vegan* 2.6.4, sobre `Data_Raw_WithCoords.csv`, com semente 123 e 9.999 permutações. Nada foi estimado, arredondado por analogia ou transportado de outro documento.

---

## 1. Scripts

### `02_constrained_ordination_and_partitioning.R`

**Três correções, cada uma documentada no ponto em que se aplica.**

**a) Frações da partição de variação lidas por posição.** Causa raiz confirmada: em *vegan* 2.6.4 as linhas de `varpart()$part$indfract` chamam-se `[a] = X1|X2`, `[b] = X2|X1`, `[c]`, `[d] = Residuals`. A linha 2 é o espaço puro e a linha 3 é a fração compartilhada. O script lia `ind[1]`, `ind[2]`, `ind[3]` e rotulava a segunda como compartilhada, o que trocou as duas em toda tabela exportada.

Verificação executada antes e depois:

```
leitura antiga  [a]=0,5162 [b]=0,0123 [c]=0,0992 | [a]+[b]=0,5285 vs adjR²(X1)=0,6154  -> FALHA
leitura nova    [a]=0,5162 compartilhada=0,0992 [c]=0,0123 [d]=0,3723
                [a]+compartilhada = 0,615383 = adjR²(X1)  -> PASSA
                compartilhada+[c] = 0,111484 = adjR²(X2)  -> PASSA
```

As frações passaram a ser lidas pelo nome de linha que o próprio *vegan* escreve, e as duas identidades aritméticas são verificadas com `stop()` antes de qualquer exportação. Se uma versão futura do *vegan* mudar a ordem, o script interrompe em vez de exportar valores trocados.

**b) Resposta da ordenação carregava o tamanho da janela.** A resposta era `sqrt` da área absoluta de classe. O primeiro componente principal dessa resposta explica 72,5% da inércia total e correlaciona-se a 0,993 com a raiz da área da unidade, que é a variável resposta do modelo univariado. Passou a ser a transformação de Hellinger, a raiz da área relativa (Legendre & Gallagher, 2001); sob ela a mesma correlação é 0,054. As duas versões são ajustadas e comparadas em `TableII_4b_RDA_Transformation_Comparison.csv`, e o diagnóstico de tamanho é exportado em `TableII_5_Response_Size_Dependence.csv`.

**c) Testes por termo eram sequenciais e reportados como marginais.** `by = "terms"` foi substituído por `by = "margin"` como teste reportado. A decomposição sequencial continua a ser exportada, agora com o nome `TableII_4_RDA_Anova_Terms_Sequential.csv`, para descrição.

**Novo:** partição sem `PD_Forest` (`TableII_3b_VP_Fractions_Without_PD.csv`), porque essa é a única métrica do conjunto paisagístico que é função algébrica da resposta.

### `03_univariate_model_selection.R`

- O bloco de entrada aceitava três arquivos em ordem de preferência, dos quais só `Data_Raw_WithCoords.csv` reproduz os resultados publicados. Agora exige esse arquivo e interrompe com instrução se ele faltar.
- Comentário corrigido: são 2⁵ = 32 subconjuntos candidatos, não 2⁴ = 16.
- Nota acrescentada acima de `cv_rmse_loo`: as funções devolvem a raiz do erro quadrático médio, e não erro absoluto médio. Foi essa confusão que entrou no manuscrito.
- `sessionInfo()` passou a ser gravado em arquivo.

### `04_geometric_coupling_diagnostic.R`

- `sessionInfo()` acrescentado e gravado em arquivo.

### `05_ratio_artefact_test.R`

- A função `emp_p()` calcula um valor de p empírico bilateral por comparação de valores absolutos, e o manuscrito o descrevia como percentil. Agora o script exporta as duas quantidades em colunas separadas, `Percentile_of_observed` e `Empirical_p_two_sided`, e a regra de decisão está escrita por extenso no cabeçalho da seção: um coeficiente é indistinguível do artefato quando cai dentro do intervalo central de 95% da nula. Uma coluna `Rules_agree` sinaliza qualquer termo em que as duas leituras discordem, que é o caso de PLAND.
- **Novo, Seção 3a:** o modelo da matriz agropecuária livre de acoplamento (`TableS_PartV_Matrix_Uncoupled.csv`). É a Tabela S40.
- **Novo, Seção 3b:** teste de desbaste para coincidência espacial entre unidades de táxons diferentes (`TableS_PartV_Spatial_Thinning.csv`).

### Arquivos de apoio

- `docs/TABLE_MAP.md`, novo: mapeia cada tabela e figura dos dois documentos ao arquivo que a produz. Era a lacuna que impedia um terceiro de localizar a origem de qualquer número.
- `CITATION.cff`: os cinco autores, afiliações, versão e palavras-chave. Restam três marcadores explícitos, listados na Seção 4 abaixo.
- `README.md`: mediana corrigida para 54,9%, Figura 7 para Figura 8, seção nova sobre o que falta no repositório e seção nova com o registro desta revisão.

---

## 2. Manuscrito

### Valores corrigidos

| Onde | Estava | Está |
|---|---|---|
| Resumo, Seção 2.3, README | mediana de 54,4% de floresta | **54,9%** |
| Resumo, Seção 3.4 | 1,83 vezes maior (complexidade de forma) | **1,82** |
| Seção 3.1 | 0,64 manchas por 100 ha em *Brachyteles* | **0,63** |
| Seção 3.1 | herbácea + não vegetada < 3% em todas as espécies | **em oito dos nove táxons; 11,2% em *L. wiedii*** |
| Seção 3.2 | "the number of patches (0.538)" | **densidade de manchas** |
| Seção 3.3, 4.1 | compartilhada 1,2%, espacial pura 9,9% | **compartilhada 9,9%, espacial pura 1,2%** |
| Seção 3.4, duas vezes | "mean absolute error" | **root mean squared error** |
| Seção 3.4 | "55.1% of the deviance", junto à Tabela 3 | **55,1% do modelo acoplado; 36,1% do modelo reportado** |
| Seção 3.4 | percentis 2,9 / 7,5 / 15,5 | **p empíricos bilaterais 0,029 / 0,075 / 0,155, com a regra de decisão declarada** |
| Seção 3.4 | "across all six specifications ... never included zero" | **cinco das seis, com a exceção nomeada** |
| Seção 2.3 | r = 0,94 atribuído à porcentagem de floresta | **r = 0,94 é da área absoluta; r = 0,04 é o da porcentagem** |
| Seção 2.4 | corpos d'água excluídos, unidades só terrestres | **mascarados os grandes corpos; água residual mediana 0,2%, máxima 10,4%** |

### Seção 3.3 reescrita

A ordenação canônica passou a ser a de Hellinger com testes marginais. A Tabela 2 foi refeita inteira, com nove linhas de dados, e a Figura 5 e a Tabela 2 tiveram as legendas atualizadas. Os valores novos, todos de execução em *vegan* 2.6.4 com 9.999 permutações:

| Preditor | Variância | F | p |
|---|---|---|---|
| Isolamento florestal (ENN_MN-F) | 0,0144 | 91,14 | 0,0001 |
| Densidade de manchas florestais (PD-F) | 0,0094 | 59,87 | 0,0001 |
| Densidade de manchas da matriz (PD-A) | 0,0027 | 17,42 | 0,0001 |
| Área média das manchas florestais (AREA_MN-F) | 0,0017 | 10,94 | 0,0006 |
| Complexidade de forma (FRAC_MN-F) | 0,0011 | 6,85 | 0,0087 |
| Proximidade da matriz (PROX_MN-A) | 0,0009 | 5,61 | 0,0135 |
| Elevação média | 0,0005 | 3,07 | 0,0756 |
| Densidade populacional humana | 0,0002 | 1,56 | 0,2201 |
| Proximidade florestal (PROX_MN-F) | 0,0001 | 0,67 | 0,4414 |

Proporção restrita 91,8%, R² ajustado 0,905, p global 0,0001. Primeiro eixo 83,6% da variância total, segundo 6,4%. Os escores do primeiro eixo correlacionam-se a 0,996 com a porcentagem de floresta.

Duas consequências foram escritas no texto sem dramatizar. A proximidade florestal continua fortemente alinhada ao gradiente composicional, com escore de biplot 0,68 no primeiro eixo, mas não carrega variação que o isolamento, a densidade de manchas e a área média já não carreguem. Isso não afeta o modelo de extensão do envelope, que é o resultado central do artigo, e a Seção 4.2 foi ajustada para dizer exatamente isso.

Retirei os valores de F por eixo do texto. Eles dependem da convenção de graus de liberdade da versão de *vegan*, e não consegui reproduzir os publicados. O texto passou a relatar a participação de cada eixo na variância total, que é independente de versão.

### Seção 4.3, a matriz agropecuária

O argumento repousava no ordenamento da RDA. Passou a repousar num modelo Gamma com ligação logarítmica que não contém nenhuma métrica normalizada pela resposta, ajustado às 67 unidades:

| Preditor | β | IC 95% | p | Razão por DP |
|---|---|---|---|---|
| Proximidade florestal | 0,893 | 0,690 a 1,095 | < 0,001 | 2,44 |
| Complexidade de forma | 0,475 | 0,278 a 0,672 | < 0,001 | 1,61 |
| **Proximidade da matriz** | **0,400** | **0,234 a 0,565** | **< 0,001** | **1,49** |
| Isolamento | 0,129 | −0,031 a 0,288 | 0,118 | 1,14 |

Deviance explicada 53,1%, VIF máximo 1,96. A densidade de manchas da matriz, no lugar da proximidade, não mostra associação (β = −0,024, p = 0,807). É a Tabela S40. O resumo passou a citar o 1,49 em vez do ordenamento da RDA.

### Outras correções de texto

- Parágrafo corrompido da Seção 2.5.1 reescrito, com a frase duplicada removida e as citações de pacotes de R acrescentadas.
- `betadisper` passou a ser descrito como PERMDISP e não PERMANOVA.
- Citação de Benjamini e Hochberg (1995) acrescentada onde a correção FDR é aplicada.
- Afirmação sobre truncamento do índice de proximidade substituída pelo fato: o raio de 28.230 m excede o diâmetro circular equivalente de 54 das 67 unidades.
- A leitura da proximidade como conectividade passou a declarar que, em posto, ela ordena as unidades quase como a área média das manchas, com rho de Spearman de 0,91.
- Título de Brooks et al. (2019) restaurado, que a substituição global de "Area of Habitat" havia corrompido.
- Hesselbarth et al. (2019) corrigido de *Ecology* para *Ecography*.
- Legendre e Gallagher (2001) e Cushman et al. (2008) acrescentados à lista, ambos verificados na fonte.
- "Brazil, 2006" para "Brasil, 2006"; acentos em Beltrão, Stăncioiu e Carmo Pônzio; Solórzano desambiguado para 2021a.
- Fahrig et al. (2026) passou a ser citado na Seção 4.1, onde é diretamente pertinente.
- Declaração de disponibilidade de dados redigida, com três marcadores explícitos a preencher.

---

## 3. Material suplementar

- Tabela S16 refeita com todas as colunas. Estava truncada em Model e k, e o corpo cita dela o AICc, o BIC e os dois erros de validação cruzada.
- Legenda da Tabela S19: intercepto, graus de liberdade, número de modelos e model averaging corrigidos.
- Legenda da Tabela S20: valores de importância corrigidos para os da própria tabela.
- Legenda da Tabela S21: coeficientes corrigidos de −0,355 e −0,021 para −0,642 e −0,542.
- Legenda da Tabela S23: "all VIF values remain below 2.3" substituída pelo que a tabela mostra, com remissão à Tabela S40.
- Legenda e título da Tabela S26: deixou de afirmar que M1 é o modelo do corpo.
- Abertura da Seção 3 e legendas das Figuras S5 e S6: deixaram de chamar o modelo retido de modelo misto.
- Stress da NMDS corrigido de 0,089 para 0,084.
- Cabeçalho da Seção 5 de 4.1 para 4.5; duas remissões à Figura 7 para Figura 8; Apêndice S1 de 4.2 para 4.6.
- Legenda da Tabela S36: o termo retirado "Areas of Habitat" substituído.
- Bartón alinhado a 2023, como na lista principal.
- **Três tabelas novas**, todas de execução verificada: S38, partição sem PD; S39, efeito da transformação da resposta sobre os testes marginais; S40, matriz agropecuária livre de acoplamento.

---

## 4. O que depende de você

**Antes de qualquer coisa: rodar `run_all.R` do início ao fim.** As correções nos scripts mudam as saídas da Parte II. Sem rodar, as tabelas de saída no disco continuam sendo as antigas.

1. **A Figura 5 ainda é a antiga.** O biplot embutido no manuscrito é o da RDA de `sqrt(CA)`. O script 02 corrigido gera o novo automaticamente em `PartII/Figure_RDA_Biplot.tiff`, mas é preciso rodar e substituir a imagem no Word. Mesma coisa para a Figura S8, o diagrama de Venn.

2. **Falta um cabeçalho "2.5" no manuscrito.** A Seção 2.4 é seguida por 2.5.1 sem que exista 2.5, e 2.5.2 também não existe. Não corrigi por edição de texto porque exige criar um parágrafo com o estilo de título certo. Inserir no Word um cabeçalho "2.5 Statistical analysis" acima de 2.5.1, e renumerar 2.5.3, 2.5.4 e 2.5.5 para 2.5.2, 2.5.3 e 2.5.4, atualizando as remissões nos dois documentos.

3. **Itálico em *Leopardus wiedii*** na Seção 3.1. A frase nova entrou em texto simples.

4. **Três marcadores na declaração de disponibilidade de dados:** `[REPOSITORY URL]`, `[DOI]` e `[VERSION]`. Os mesmos três estão no `CITATION.cff` junto com `date-released`. Preencha os quatro de uma vez.

5. **Duas afiliações e um telefone** continuam como `[INSERT]`.

6. **Yue et al. (2015).** A citação estava na Seção 4.3 e não tinha entrada na lista; procurei na literatura de land sparing e land sharing sem correspondência inequívoca. Retirei a citação e mantive as outras três. Se havia um trabalho pretendido, reponha com os dados completos.

7. **Rubia Morini** consta dos agradecimentos e recebe atribuição de contribuição de autoria como responsável pela visualização. Ou é autora, ou a contribuição vai para outra pessoa.

8. **Área do estado:** o manuscrito diz 43.750 km² e a auditoria interna anterior registra 43.700 km². Confira contra a fonte do IBGE e use um valor só.

9. **Coleção do MapBiomas:** os valores de acurácia, 91,5% no nível 1 e 86,1% no nível 2 para a Mata Atlântica, e a unidade mínima de mapeamento de seis pixels, estão corretos no apêndice oficial. A documentação designa a coleção como "Collection 10" e o manuscrito diz "Collection 10.1". Confirme a designação.

10. **Periódico.** O corpo tem cerca de 9.500 palavras, o resumo 301 e há 11 itens gráficos. Isso decide quanto migra para o suplementar. Se for PECON, o corte é grande. Ainda vinte e poucas referências da lista não são citadas em lugar nenhum: a limpeza depende do corte.

11. **O repositório continua sem dados.** `data/raw/` e `data/processed/` têm só `.gitkeep`. Sem eles, `run_all.R` para na primeira linha. O `README.md` agora lista os quatro conjuntos que faltam.

12. **Resposta armazenada como inteiro.** `AOH_unit_area_ha` é do tipo inteiro e 28 dos 67 valores são múltiplos de 100. Exporte do ArcGIS com precisão total e rode de novo. O efeito nos resultados deve ser desprezível, mas a precisão da resposta faz parte do argumento de um artigo sobre identidade algébrica entre resposta e preditores.

13. **Resultado do teste de desbaste.** Rodando o script 05 você obterá `TableS_PartV_Spatial_Thinning.csv`. Na minha execução, com 56 unidades por réplica e 200 réplicas, a proximidade tem mediana 0,806 e a complexidade de forma 0,658, ambas com p abaixo de 0,05 em 100% das réplicas, contra 0,741 e 0,601 no desenho completo. Ou seja, as duas associações reportadas ficam mais fortes quando a coincidência espacial é removida, e não mais fracas. Vale entrar no suplementar como tabela e ganhar uma frase na Seção 4.6, depois que você rodar e confirmar os números no seu ambiente.

---

## 5. O que não foi alterado, e por quê

- **Os testes por eixo da RDA.** Não reproduzi os valores publicados de F por eixo em *vegan* 2.6.4; a convenção de graus de liberdade difere entre versões. Em vez de substituir números que não pude confirmar, retirei os F por eixo do texto e deixei a participação de cada eixo na variância total, que não depende de versão.
- **As Tabelas S36 e S37 e a Figura 8.** Continuam sem script. A caixa de proveniência do suplementar declara isso com honestidade, e mantive.
- **A estrutura do modelo univariado.** Nada nele estava errado. A Parte V é a melhor parte do trabalho e foi reproduzida coeficiente a coeficiente, até a quarta casa decimal, fora do seu ambiente de R.
- **A Tabela 3.** Não mudou um valor. As duas associações reportadas resistiram a tudo que testei.
