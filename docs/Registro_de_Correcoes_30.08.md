# Registro de correcoes, 30 de agosto de 2026

A conferencia desta data comecou por uma discrepancia de horario e terminou
resolvendo o problema mais serio que a cadeia tinha. As saidas dos scripts 02 a
05 continuavam datadas de 25 de agosto, embora a autora tivesse rodado o
`run_all.R` no dia 29. A cadeia havia rodado; escrevera no lugar errado.

Arquivos entregues:

- `Manuscript_REVISED_30.08.docx`
- `Supplementary_Material_REVISED_30.08.docx`
- `Repositorio_GitHub/run_all.R`, com verificacao de ancora e log em arquivo

As versoes de 29 de agosto estao em `_Arquivo/Manuscritos_superados/`.

---

## 1. Por que as saidas nao apareceram

O pacote `here` procura uma ancora para definir a raiz do projeto, e o arquivo
`.here` guardado dentro de `Repositorio_GitHub/` e uma dessas ancoras. Ele existe
para que alguem que clone apenas o repositorio consiga rodar a cadeia. Contudo,
quando a sessao de R comeca dentro dessa pasta, e foi o que ocorreu em 29 de
agosto, `here()` resolve para `Repositorio_GitHub` e nao para a pasta que contem
`ANALISES_TCC.Rproj`. Nada da erro. Todo `here("Outputs", "Manuscrito")` aponta um
nivel abaixo, e os arquivos vao para `Repositorio_GitHub/Outputs/Manuscrito/`.

A prova esta nos proprios arquivos: os `session_info_part02` a `part05` foram
gravados naquela pasta em 29 de agosto, entre 18h20 e 18h22. A cadeia completou.

## 2. A comparacao entre as duas corridas

Foram comparados os 54 CSV presentes nas duas arvores, valor a valor.

| Resultado | Arquivos |
|---|---|
| Identicos em todos os valores | 51 |
| Diferentes apenas no rotulo de uma coluna | 3 |
| Diferentes em algum valor | 0 |

Os tres sao `TableS_PartV_Uncoupled_Model.csv` e
`TableS_PartV_Proximity_Conditioned.csv`, em que `AOH_ratio_per_SD` passou a
`OHE_ratio_per_SD`, e `Table_S_PartIII_DescriptiveStatistics.csv`, em que
`log_AOH` passou a `log_OHE`. Sao as renomeacoes registradas em 26 de agosto nos
scripts 03 e 05, que so agora foram executadas. O suplemento ja imprimia os
rotulos novos, de modo que a corrida de 29 de agosto e a primeira em que os dois
concordam.

Nenhum numero do manuscrito muda por causa disso. A corrida de 29 de agosto
reproduziu a de 25 de agosto em todos os valores, o que e a evidencia de
reprodutibilidade mais forte que o trabalho tem ate agora.

## 3. O que foi feito com os arquivos

As saidas de 29 de agosto foram copiadas para `Outputs/Manuscrito/`, que passa a
ser a arvore corrente. A arvore que ficara dentro do repositorio foi movida para
`_Arquivo/Outputs_gravados_no_repo_29.08/`, e `Repositorio_GitHub/Outputs/` voltou
a conter apenas o `.gitkeep`. Os tres arquivos de 25 de agosto cujos rotulos
mudaram estao em `_Arquivo/Outputs_25.08_substituidos/`.

O script 01 continua sendo a excecao. Ele gravou 13 dos seus arquivos entre 18h16
e 18h18 e parou antes da Secao 10, de modo que o `session_info_part01.txt`, as
Figuras 1 a 3 e as Tabelas 3, 3b e 4 permanecem com data de 25 de agosto. Os 13
arquivos regravados reproduziram os valores anteriores exatamente, incluindo o
estresse do NMDS de 0,0844, a PERMANOVA de genero com R2 de 0,227 e F de 2,4745, e
os sete R2 do envfit.

## 4. A protecao acrescentada ao run_all.R

O `run_all.R` agora para com mensagem explicita se `here()` resolver para uma pasta
que nao contenha `ANALISES_TCC.Rproj`, e grava `Outputs/Manuscrito/run_all_log.txt`
com data e hora de cada etapa. Cada script roda dentro de um `tryCatch` que
registra a mensagem de erro em vez de derrubar a cadeia em silencio, e o log conta
quantos arquivos cada script gravou, avisando quando esse numero for zero. As
variaveis internas do executor receberam prefixo de ponto para que nenhum script
carregado por `source()` as sobrescreva.

## 5. Edicoes no manuscrito

### 1. Figura 6, trocada

> **Antes:** o painel de seis termos do script 03, cujos valores vem da media de
> modelos da Tabela S18b (proximidade 1,363; forma 0,385; isolamento -0,352)
>
> **Depois:** `Figure_PV_02_Reported_Model_Coefficients.tiff`, que plota os tres
> coeficientes da Tabela 3 (0,742; 0,602; 0,122)
>
> Motivo: a ultima frase da legenda afirmava "Table 3 shows the model from the
> upper panel alone; Figure 6 shows its coefficients". Isso nao era verdade. O
> painel superior trazia os coeficientes medios do conjunto candidato, e a Tabela 3
> traz o ajuste unico da especificacao S4, que sao numeros diferentes. O leitor via
> 1,36 na figura e 0,742 na tabela.

### 2. Legenda da Figura 6, reescrita

> **Depois:** "Coefficients of the model reported in Table 3, a Gamma generalized
> linear model with a logarithmic link fitted to the 67 sampling units, from which
> every metric normalized by the area of the sampling window was removed. [...]
> A proximity coefficient of 0.742 means the envelope is 2.10 times larger for each
> standard deviation of clustering among remnants. Figure S10 shows the full a
> priori set of six terms, three of which the permutation null reproduces and which
> the paper therefore does not interpret."

### 3. Resolucao das figuras

As oito figuras estavam embutidas entre 150 e 220 dpi, com 973 a 1.473 pixels de
largura, enquanto as do suplemento estavam a 600 dpi. A Figura 2 ainda carregava
uma etiqueta de resolucao de 1 dpi. As Figuras 2 a 7 foram substituidas pelos TIFF
exportados, a 600 dpi, e a extensao da Figura 6 na pagina foi ajustada porque a
proporcao mudou de 1,612 para 1,751 com a troca.

As Figuras 1 e 8 continuam a 150 dpi. Sao mapas do ArcGIS e nao existe versao de
alta resolucao em lugar nenhum da pasta. Precisam ser reexportadas a 600 dpi
antes da submissao.

### 4. Secao 4.1, frase repetida palavra por palavra

> **Antes:** Model selection is routinely treated as a neutral arbiter. However, here it is not: a predictor built from the response improves the fit to the response whether or not it carries ecological information, so the criterion measures access to the answer rather than explanatory value.
>
> **Depois:** Model selection is routinely treated as a neutral arbiter, and Section 3.4 shows that here it is not, because the criterion measures access to the answer rather than explanatory value.
>
> Motivo: a oracao removida aparecia identica na Secao 3.4. A Discussao interpreta
> o resultado; nao precisa reapresenta-lo.

### 5. Secao 4.1, fracoes da particao repetidas

> **Antes:** ... (51.6%), and 47.5% once the one metric normalized by the response is removed from the set. Two Moran eigenvectors were retained, indicating spatial structure. With the full landscape set, the purely spatial fraction is 1.3% and not significant, while 9.9% is shared with the landscape. However, without that metric, the purely spatial fraction rises to 6.8% and becomes significant ...
>
> **Depois:** Landscape structure nonetheless accounted for a substantial share of the variation in habitat amount that is independent of geography, and Section 3.3 reports how that share changes once the one metric normalized by the response is removed from the set.
>
> Motivo: as seis porcentagens ja estao na Secao 3.3, na mesma ordem e com os
> mesmos testes. O paragrafo perde 40 palavras e mantem o que so ele diz, que sao
> os 37,2% nao explicados e o controle espacial da Secao 3.4.

### 6. Secao 4.6, terceira ocorrencia da mesma exigencia

> **Antes:** Separating composition from configuration in the sense of Fahrig (2013) requires landscape metrics measured in windows defined independently of the species data, on a fixed grid or in buffers of fixed radius around occurrence points, as in Rios et al. (2021), and that is the priority for further work in this region.
>
> **Depois:** Separating composition from configuration in the sense of Fahrig (2013) requires the design set out in Section 4.1, and that is the priority for further work in this region.
>
> Motivo: a exigencia e a citacao de Rios et al. (2021) ja aparecem na Secao 4.1 e
> na Conclusao. Tres enunciados da mesma condicao em um artigo do tamanho deste
> cansam o leitor sem acrescentar evidencia.

O corpo passou de 10.874 para 10.619 palavras.

### 6b. Resumo, precisao restaurada

O resumo havia sido reescrito pela autora entre as sessoes, passando de 297 para
185 palavras. A reducao melhorou o ritmo, mas retirou cinco elementos e enfraqueceu
uma afirmacao.

> **Antes:** Two associations remained significant: larger envelopes correlated with clustering (2.10x) and complex patch shapes (1.83x), even after accounting for occurrence records. [...] the null model produced similar results, suggesting no ecological effect.
>
> **Depois:** Two associations survived it: per standard deviation, envelopes were 2.10 times larger with clustering among remnants and 1.83 times larger with patch shape complexity, falling to 1.51 and 1.63 times after controlling for the number of occurrence records inside the envelope. [...] the null reproduced it, and also the coefficients of forest cover and patch density, so none of the three is distinguishable from the value the ratio construction generates on its own.
>
> Motivo: uma razao sem a unidade a que se refere nao informa nada, e 2,10 so tem
> sentido por desvio padrao. As razoes condicionais de 1,51 e 1,63 precisam aparecer
> ao lado das incondicionais, ou o leitor toma 2,10 e 1,83 como se nao dependessem
> do esforco amostral. E dizer que a nula sugere "no ecological effect" afirma mais
> do que o artigo sustenta: a nula torna o coeficiente indistinguivel do artefato de
> razao, o que nao e o mesmo que demonstrar ausencia de efeito.

Foram tambem devolvidos os 78,8% da mancha unica de *Brachyteles*, o numero de
nove taxons, e o italico dos nomes de genero. O resumo ficou com 260 palavras.

## 6. Edicao no material suplementar

### 7. Figura S10, nova

O painel de seis termos que saiu do corpo entrou no suplemento como Figura S10, a
600 dpi, com legenda que declara de onde vem cada metade do painel e por que a
metade inferior e desenhada sem marcas de significancia. A numeracao das Figuras
S1 a S9 nao mudou.

## 7. Conferencia numerica

Manuscrito: 3 tabelas e 11 legendas, nenhum numero sem contrapartida em CSV.
Figuras 1 a 8 e Tabelas 1 a 3 em sequencia completa.

Suplemento: 55 tabelas e 60 legendas. Tabelas S1 a S42 e Figuras S1 a S10 em
sequencia completa. Quatro tabelas contem numeros sem contrapartida em CSV, e as
quatro sao declaradas no proprio texto: o raio de busca de 28.230 m na Tabela S2,
os 606 registros retidos na Tabela S42b, e as Tabelas S36 e S37 do indice de valor
de conservacao, que vem do ArcGIS.

## 8. O que continua pendente

1. Reexportar as Figuras 1 e 8 do ArcGIS a 600 dpi.
2. Rodar o script 01 ate o fim, para que o `session_info_part01.txt` e as Figuras
   1 a 3 do PartI passem a vir da mesma corrida que o restante da cadeia.
3. `git init` no `Repositorio_GitHub/`, que ainda nao e um repositorio git.
4. Linhas 394 a 400 do script 03, que aceitam tres arquivos de entrada onde so um
   produz os resultados publicados.
5. Depositar os 67 rasters, 2,6 GB, em repositorio com DOI e registrar o DOI no
   `CITATION.cff`, hoje com tres marcadores vazios.
6. Duas afiliacoes, telefone do autor correspondente, declaracao de
   disponibilidade de dados.
7. A atribuicao de Rubia Morini, hoje nos agradecimentos e na contribuicao de
   autoria ao mesmo tempo.
