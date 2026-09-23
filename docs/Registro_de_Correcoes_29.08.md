# Registro de correcoes, 29 de agosto de 2026

Trabalho feito depois da execucao do script 08 no R e de uma auditoria numerica
por rastreamento, nao por correspondencia. Cada valor impresso no manuscrito e no
suplemento foi procurado nos 69 arquivos CSV exportados pela cadeia, em todos os
arredondamentos plausiveis, e cada celula de cada uma das 58 tabelas dos dois
documentos foi submetida ao mesmo teste.

Arquivos entregues:

- `Manuscript_REVISED_29.08.docx`
- `Supplementary_Material_REVISED_29.08.docx`
- `Guia_Reuniao_Orientador_29-08-2026.docx`
- `Repositorio_GitHub/Dados/Processados/Data_Raw_WithCoords.csv`, substituido

As versoes de 26 de agosto estao em `_Arquivo/Manuscritos_superados/`.

---

## 1. O que a auditoria encontrou

O resultado foi limpo. Das 58 tabelas dos dois documentos, apenas quatro contem
numeros sem contrapartida em CSV, e as quatro sao declaradas: a Tabela S2, que
imprime o raio de busca de 28.230 m; a Tabela S42b, que traz os 606 registros
retidos antes do recorte por unidade; e as Tabelas S36 e S37, do indice de valor
de conservacao, que vem do ArcGIS e nao da cadeia em R, como a propria caixa de
proveniencia do suplemento declara.

No texto corrido, todos os numeros sinalizados eram intervalos de pagina da lista
de referencias, exceto quatro constantes externas: a area do estado, a area de
referencia de *Panthera onca*, o raio de busca e a extensao mapeada do indice de
valor de conservacao.

Uma unica discrepancia real apareceu, e foi corrigida.

---

## 2. Edicoes no manuscrito

### 1. Secao 3.4, deviance explicada pelo modelo reportado

> **Antes:** ... accounted for 36.1%.
>
> **Depois:** ... accounted for 36.0%.
>
> Motivo: a unica fonte exportada para essa quantidade e
> `Outputs/Manuscrito/PartVI/TableS42_Sampling_Effort.csv`, que registra 0,3604
> para o modelo M1, o modelo da Tabela 3. O valor de 0,360784 citado no registro
> de 26 de agosto veio da implementacao independente em Python e nao tem
> contrapartida em nenhum CSV atual. Os tres scripts que calculam a quantidade
> usam a mesma formula, `1 - deviance / null.deviance`, de modo que a diferenca
> esta na implementacao e nao no metodo. Com a correcao, todo numero do corpo do
> manuscrito passa a ter origem rastreavel num arquivo exportado.

### 2. Resumo, de 340 para 297 palavras

Nenhum dos doze numeros do resumo foi perdido: 67 unidades, mediana de 54,9%,
18 de 67 unidades abaixo de 30%, *Mazama* a 28,5%, *Alouatta* a 29,3%,
*Brachyteles* a 81,2% com mancha unica de 78,8%, as razoes de 2,10 e 1,83, as
razoes condicionais de 1,51 e 1,63, e a razao de 1,49 da matriz consolidada.

As reducoes foram de redundancia e nao de conteudo. Tres exemplos.

> **Antes:** with notable variation: 18 of the 67 units contained less than 30% forest, a threshold linked to abrupt declines in this biome
>
> **Depois:** and 18 of the 67 units held less than 30%, below which biodiversity declines abruptly in this biome
>
> Motivo: "with notable variation" nao acrescenta informacao a uma frase que ja
> exibe a variacao, e "a threshold linked to" e mais longo e mais vago que a
> oracao relativa.

> **Antes:** the envelope was 2.10 times larger per standard deviation increase in clustering among forest remnants, and 1.83 times larger per standard deviation increase in patch shape complexity
>
> **Depois:** per standard deviation, the envelope was 2.10 times larger with clustering among forest remnants and 1.83 times larger with patch shape complexity
>
> Motivo: "per standard deviation" aparecia duas vezes na mesma frase.

> **Antes:** Envelopes were also larger where the farming matrix itself was spatially consolidated rather than dispersed, by a factor of 1.49 per standard deviation
>
> **Depois:** Envelopes were also 1.49 times larger where the farming matrix was spatially consolidated rather than dispersed
>
> Motivo: "by a factor of" repete o que "times larger" ja diz, e "itself" nao
> tem antecedente contrastivo na frase.

### 3. Secao 4.4, densidade populacional humana, de 111 para 100 palavras

> **Antes:** Human population density predicts neither composition nor habitat amount, which is more likely a property of the metric than an absence of pressure. A density averaged across an entire sampling unit summarizes residence rather than activity and obscures local gradients of pressure on these populations.
>
> **Depois:** Human population density predicted neither composition nor habitat amount, which is more likely a property of the metric than an absence of pressure: a density averaged over an entire sampling unit summarizes residence rather than activity.
>
> Motivo: a secao relata um resultado nulo de uma covariavel que nao sustenta
> nenhuma conclusao. As duas frases diziam a mesma coisa, e "obscures local
> gradients of pressure" e outra formulacao de "summarizes residence rather than
> activity". O tempo verbal tambem foi alinhado ao restante dos Resultados.

---

## 3. Edicao no material suplementar

### 4. Cabecalho duplicado na Secao 1

> **Antes:** dois cabecalhos numerados 1.2, um deles "1.2 Landscape Metrics: Definitions and Ecological Interpretation" e o outro "1.2 Landscape Metrics", este ultimo servindo apenas de guarda-chuva para "1.2.1 Technical Notes on Metric Calculation"
>
> **Depois:** o segundo cabecalho foi removido, e "1.2.1 Technical Notes on Metric Calculation" passa a ser subsecao do primeiro
>
> Motivo: a numeracao duplicada estava registrada como pendencia desde o parecer
> de 22 de agosto. Remover o cabecalho redundante resolve a duplicacao sem
> renumerar nada a jusante.

---

## 4. Correcao no repositorio

### 5. `Repositorio_GitHub/Dados/Processados/Data_Raw_WithCoords.csv`

A copia do repositorio trazia a variavel resposta arredondada para inteiro, com
valores diferentes dos publicados. A primeira unidade de *Alouatta guariba*
constava com 126.805 ha onde a versao corrente tem 126.107,31 ha, e a segunda com
16.834 ha onde a versao corrente tem 16.765,21 ha. Quem clonasse o repositorio
obteria numeros diferentes dos do artigo, sem aviso.

O arquivo foi substituido pela versao de `Dados/Processados/`. Os dois md5 agora
conferem: `881996be7335b4af0c90bcc5c5303e41`.

---

## 5. O que continua pendente

1. A pasta `Repositorio_GitHub/` ainda nao e um repositorio git. Sem `git init` e
   um primeiro commit nao ha o que publicar no GitHub.
2. O script 03, nas linhas 394 a 400, ainda aceita tres arquivos de entrada
   alternativos, dos quais so `Data_Raw_WithCoords.csv` produz os resultados
   publicados.
3. Depositar os 67 rasters por unidade, 2,6 GB, num repositorio com DOI e
   registrar o DOI no `CITATION.cff`.
4. Preencher `date-released`, `repository-code` e `doi` no `CITATION.cff`.
5. Duas afiliacoes, o telefone do autor correspondente e a declaracao de
   disponibilidade de dados.
6. A atribuicao de Rubia Morini, que consta dos agradecimentos e recebe
   atribuicao de contribuicao de autoria.

## 6. Redundancias identificadas e nao removidas

Tres repeticoes de conteudo sobrevivem no corpo, e a decisao sobre elas foi
adiada. Ficam registradas com a localizacao exata.

1. A frase "a predictor built from the response improves the fit to the response
   whether or not it carries ecological information" aparece na Secao 3.4 e de
   novo na Secao 4.1, praticamente palavra por palavra. A ocorrencia da Secao 4.1
   pode sair sem perda.
2. As fracoes da particao de variacao, 51,6%, 47,5% e 37,2%, sao reportadas na
   Secao 3.3 e repetidas na Secao 4.1 com redacao quase identica. A Secao 4.1
   poderia remeter a Secao 3.3 em vez de reapresenta-las.
3. A exigencia de "habitat amount measured over windows defined independently of
   the species data" aparece na Secao 4.1 e na Conclusao. Uma das duas basta.

Removendo as tres, o corpo perde cerca de 130 palavras e nenhum argumento.
