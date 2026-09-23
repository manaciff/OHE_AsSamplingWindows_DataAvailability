# Registro de correções, 25 de agosto de 2026

Trabalho feito depois da reexecução completa da cadeia pela autora e da sua
segunda passagem de linguagem sobre o manuscrito e o material suplementar.

Arquivos entregues:

- `Manuscript_REVISED_25.08.docx`
- `Supplementary_Material_REVISED_25.08.docx`
- `Figura_Conectividade_dPC.png` e `TabelaS_Conectividade_dPC_demonstracao.csv`
  (demonstração para a defesa, fora do manuscrito)
- Guia interativo de defesa, publicado como artefato

---

## 1. O que foi conferido

A cadeia reexecutada gravou 66 arquivos CSV em `Outputs/Manuscrito/`. Comparados
um a um com a execução anterior, apenas três diferem:

| Arquivo | O que mudou |
|---|---|
| `TableS28_Geometric_Coupling.csv` | todos os valores, na quarta casa decimal |
| `PartV/TableS_PartV_Identity_and_PartialCorrelations.csv` | apenas rótulos, AOH → OHE; os valores são idênticos |
| `Table_S_PartIII_ObsPred_R2.csv` | arquivo novo, exportado pela correção de 24 de agosto |

Os outros 63 são idênticos byte a byte, de modo que a conferência anterior
continua valendo para eles.

### A pendência de 24 de agosto está resolvida

O script 04 lia `Data_Raw_FINAL.csv` enquanto os scripts 03 e 05 liam
`Data_Raw_WithCoords.csv`. A correção foi aplicada em 24 de agosto e a
reexecução confirma o efeito: a verificação da identidade algébrica passou de
**0,029 para 0,0277**, que é exatamente o valor que o script 05 calcula de forma
independente. A Tabela S26 do suplemento já trazia os valores corretos e agora
concorda com o arquivo que os produz.

### Verificação célula a célula

Foi escrito `doc/verify_grid2.py`, que para cada tabela do suplemento:

1. casa cada linha do documento com a linha correspondente do CSV pela chave
   (uma ou duas colunas), de modo que um valor não pode passar por estar em
   outro lugar do arquivo;
2. para cada coluna numérica, exige que exista pelo menos um campo do CSV que
   reproduza o valor impresso **em todas as linhas**, na precisão impressa.

Resultado: **38 de 42 tabelas aprovadas automaticamente**. As quatro restantes
foram conferidas à mão e também estão corretas:

- **S33**: chave composta com formatação diferente; conferida por script próprio,
  26 linhas, nenhuma divergência.
- **S38**: a coluna "% do total" é derivada (R² ajustado × 100) e por isso não
  corresponde a um campo único do CSV.
- **S39**: as duas metades da tabela vêm de dois grupos do mesmo arquivo
  (`Response` = sqrt(class area) e Hellinger); os 36 valores foram conferidos um
  a um.
- **S40**: a linha do intercepto não é exportada pela cadeia. O modelo foi
  reajustado em R sobre `Data_Raw_WithCoords.csv` e reproduziu
  10,5603 / 0,0732 / [10,4168; 10,7037], além de 53,07% de deviance explicada e
  VIF máximo 1,963.

As três tabelas do manuscrito (Tabelas 1, 2 e 3) foram conferidas contra
`TableS4_PERMANOVA_Global.csv`, `TableII_4_RDA_Anova_Marginal.csv` e
`TableS_PartV_Uncoupled_Model.csv`. Nenhuma divergência.

---

## 2. O que foi corrigido

### 2.1 Manuscrito, Seção 3.4: duas frases que contradiziam a tabela citada

A segunda passagem de linguagem reintroduziu duas inversões que já haviam sido
corrigidas em 24 de agosto.

> **Antes:** "The clustering coefficient ranged from 0.67 to 1.43 and **was
> excluded from** five of the six specifications."
>
> **Depois:** "…and **excluded zero in** five of the six specifications."
>
> Motivo: a Tabela S28a mostra proximidade presente nas seis especificações. O
> que ocorre em cinco delas é o intervalo de confiança excluir o zero. A redação
> anterior dizia que o preditor fora removido do modelo.

> **Antes:** "…the absolute class area of forest, **a quantity that bounds the
> response from above**…"
>
> **Depois:** "…**a quantity the response bounds from above**…"
>
> Motivo: a floresta dentro de um envelope não pode ser maior que o envelope,
> então a resposta limita a área de classe, e não o contrário. A Tabela S28e diz
> isso na direção correta, de modo que os dois documentos discordavam.

### 2.2 Suplementar, legenda da Tabela S21: dois arredondamentos errados

> **Antes:** "…is 0.659 (Pearson) and 0.723 (Spearman)."
>
> **Depois:** "…is 0.660 (Pearson) and 0.729 (Spearman)."
>
> Motivo: `Table_S_PartIII_StabilityCheck.csv` traz 0,6598 e 0,7289, impressos na
> própria tabela com quatro casas. A legenda arredondava mal a primeira e errava
> a terceira casa da segunda.

### 2.3 Suplementar, legenda da Tabela S28f: oração sem predicado

> **Antes:** "…the mean forest patch area, the one measure of amount here that
> **the response has no bounds**."
>
> **Depois:** "…that **the response neither contains nor bounds**."
>
> Motivo: a passagem de linguagem removeu o verbo e, com ele, o sentido.

### 2.4 Correções tipográficas e de formatação

- "Section 4. 4.6." → "Section 4.6." (Apêndice S1)
- "a 30- year window" → "a 30-year window" (Apêndice S1)
- Itálico restaurado em quatro trechos em que o Word uniu execuções e perdeu a
  marcação: três no manuscrito, um no suplementar (§5 do suplementar, parágrafo
  dos pesos de sensibilidade).

### 2.5 O que **não** foi alterado

- Os cinco marcadores amarelos das afiliações no manuscrito, por decisão da
  autora.
- Nenhum valor numérico de nenhuma tabela, porque nenhum estava errado.
- A Figura S10: a autora substituiu o arquivo pelo
  `Figure_PIII_02_Coefficient_ForestPlot.tiff` recém-gerado pelo script 03.
  Verificado: 4110 × 2551 px a 600 dpi, com os dois painéis, os três termos
  acoplados em cinza e ocos, e sem marca de significância. A legenda existente
  descreve exatamente essa figura.

---

## 3. Demonstração do índice de conectividade

Fora do manuscrito, a pedido da autora e depois de consulta.

### 3.1 Por que é aplicável

A probabilidade de conectividade de Saura e Pascual-Hortal (2007) é

    PC = (1 / A_L²) · Σ_i Σ_j a_i · a_j · p*_ij

com `A_L` igual à área da paisagem, que neste desenho é a área do envelope, ou
seja, a variável resposta. **PC herda o acoplamento da Seção 2.5.1 e não pode ser
comparado entre unidades amostrais.**

A importância de uma mancha, porém, é uma razão:

    dPC_k = 100 · (PC − PC_sem_k) / PC

e `A_L²` aparece no numerador e no denominador. Remover uma mancha não altera
`A_L`, então o termo cancela exatamente. **O ranking de manchas dentro de um
mesmo envelope é livre do acoplamento.** É essa distinção que torna o mapa de
prioridades defensável e uma regressão sobre PC não.

### 3.2 O que foi calculado

Três unidades, escolhidas para contrastar. Manchas de floresta extraídas do
raster de 30 m com regra de oito vizinhos, distâncias borda a borda entre pares,
`p_ij = exp(−k·d_ij)` calibrado para 0,5 na distância mediana de dispersão,
ligações cortadas acima de três vezes essa distância, e `p*_ij` obtido pelo
caminho de produto máximo. As três distâncias medianas, 250, 1.000 e 2.500 m,
formam uma escada de sensibilidade, porque a distância de dispersão **não é
medida aqui**.

| Unidade | Manchas | Floresta (ha) | Envelope (ha) | Manchas com metade do dPC | Soma da fração conector | rho(dPC, área) |
|---|---:|---:|---:|---:|---:|---:|
| *Brachyteles arachnoides* UA 1 | 39 | 22.388 | 24.613 | 1 | 0,0% | 0,987 |
| *Leopardus wiedii* UA 1 | 175 | 86.488 | 95.121 | 1 | 0,0% | 0,990 |
| *Mazama* UA 1 | 1.082 | 13.702 | 50.669 | 7 | 70,8% | 0,958 |

Valores com distância mediana de 1.000 m.

### 3.3 O que a demonstração mostra

Nos dois envelopes ricos em floresta o resultado é trivial: uma única mancha
responde por praticamente todo o índice, e nenhuma outra funciona como
trampolim. Em *Mazama* a prioridade se espalha, e a decomposição de Saura e
Rubio (2010) identifica trampolins. Com dispersão mediana de 250 m, a mancha 191,
de **42,5 ha**, ocupa o **sexto lugar** em dPC total (14,17), com fração conector
de 13,14 contra 0,005 de fração intra: ela quase não contribui como habitat e
contribui como conexão.

Esse é o apoio empírico que a segunda recomendação da Seção 4.5 do manuscrito
hoje não tem, e que sustenta em citações.

A dependência da suposição é forte e vale declarar: em *Mazama*, a soma da fração
conector cai de 96,3% a 250 m para 70,8% a 1.000 m e 37,0% a 2.500 m, enquanto a
correlação de posto entre dPC e área sobe de 0,753 para 0,991. Quanto maior a
distância de dispersão assumida, mais o índice vira uma medida de tamanho de
mancha.

### 3.4 Por que não entra no manuscrito

Três razões, e a decisão foi da autora.

1. Exigiria distâncias de dispersão por táxon que estes dados não têm. A relação
   alométrica de Sutherland et al. (2000) daria 4,9 km para *Brachyteles* e
   7,3 km para *Mazama*, valores em que o índice já colapsou para área.
2. Responde a outra pergunta: quais manchas priorizar **dentro** de um envelope,
   e não quais atributos da floresta acompanham envelopes maiores.
3. Acrescentar uma análise nova a poucos dias da defesa é risco sem retorno
   proporcional.

### 3.5 Arquivos

- `R/06_connectivity_demo.R` — funções de extração de manchas, distâncias e dPC
- `R/run_demo3.R` — executa as três unidades
- `R/07_connectivity_map.R` — mapa e tabela
- `Manuscrito/Guias_DocsSuplementares_Estudos/Figura_Conectividade_dPC.png`
- `Manuscrito/Guias_DocsSuplementares_Estudos/TabelaS_Conectividade_dPC_demonstracao.csv`

Nenhum desses scripts é chamado por `run_all.R`. A cadeia inferencial permanece
com seis scripts.

---

## 4. Duas inconsistências de rótulo nos arquivos exportados

Não afetam nenhum valor, e ficam registradas para a próxima passagem nos scripts:

- `Table_S_PartIII_DescriptiveStatistics.csv` ainda escreve `AOH_ha` e `log_AOH`,
  enquanto a Tabela S15 do suplemento imprime `OHE_ha` e `log_OHE`.
- `PartV/TableS_PartV_Uncoupled_Model.csv` ainda escreve a coluna
  `AOH_ratio_per_SD`, enquanto `TableS_PartV_Matrix_Uncoupled.csv`, do mesmo
  script, já escreve `OHE_ratio_per_SD`.

---

## 5. O que continua pendente

1. Depositar os 67 rasters por unidade, 2,6 GB, em repositório com DOI, e
   registrar o DOI no `CITATION.cff`.
2. Preencher `date-released`, `repository-code` e `doi` no `CITATION.cff`.
3. Decidir as afiliações dos coautores 2 e 3, hoje marcadas em amarelo.
