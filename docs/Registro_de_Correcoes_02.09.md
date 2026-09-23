# Registro de correções — 02 de setembro de 2026

Conferência final. O objetivo foi rastrear cada número do corpo e do material
suplementar até a célula do arquivo exportado que o produz, revisar as figuras e
as legendas, e remover o que era redundante. Os scripts já haviam sido
executados: `Outputs/Manuscrito/run_all_log.txt` registra a cadeia completa
(02, 03, 04, 05 e 08) em 02/09, das 10:10:27 às 10:12:31, sem falha.

Nenhum valor citado aqui foi estimado. Cada um foi lido do CSV indicado.

---

## 1. O que confere

As Tabelas 1, 2 e 3 do corpo e as Tabelas S3 a S11, S15 a S24, S26, S28a a S28f,
S29 a S35 e S38 a S41 do suplementar batem com os arquivos exportados hoje.

Refiz a RDA por fora da cadeia, em Python, a partir de
`Dados/Processados/Data_Raw_WithCoords.csv`, e reproduzi os valores publicados:
inércia restrita 0,9178, inércias dos três eixos canônicos (0,091364, 0,006974 e
0,002006), r = 0,996 entre os escores dos sítios e PLAND\_Forest, e escore de
biplot 0,68 da proximidade florestal no primeiro eixo. Os dois últimos não têm
CSV correspondente e eram, até aqui, os únicos números do corpo sem
rastreabilidade direta.

Também recalculei a partir dos dados brutos, e conferem: mediana de 54,85% de
cobertura florestal, 18 unidades abaixo de 30%, água residual com mediana de
0,22% e máximo de 10,36%, LPI mediana de 78,79% em *Brachyteles*, r = 0,939
entre área de classe florestal e área da unidade contra r = 0,035 para o
percentual, e 54 das 67 unidades com diâmetro circular equivalente menor que o
raio de busca de 28.230 m. As contagens das PERMANOVA par a par conferem: 3 de
28 pares de gêneros significativos após correção FDR, e 9 de 36 pares de
espécies significativos antes da correção e nenhum depois.

## 2. Três erros que contradiziam o próprio argumento do artigo

**Legenda da Figura 6.** Dizia "with all metrics normalized by area". É o
contrário: essa é a especificação que não contém métrica alguma normalizada pela
área da janela, e é sobre isso que a Seção 3.4 se apoia. Passou a "containing no
metric normalized by the area of the sampling window".

**Legenda da Figura S5, o valor.** Dizia "0.662 on the hectare scale and 0.662 on
the logarithmic scale, differing due to..." — dois valores iguais declarados
diferentes. `Table_S_PartIII_ObsPred_R2.csv` traz 0,6623 na escala de hectares e
0,5570 na escala logarítmica. O segundo passou a **0,557**.

**Legenda da Figura S5, o modelo.** Atribuía o painel ao "Gamma model ... as
shown in Table 3". O script 03 constrói o painel a partir de `predict_model`, que
é o modelo acoplado; o próprio CSV registra `Model used for the figure: Gamma,
interaction`. Atribuí-lo à Tabela 3 contradizia a Tabela S35, que dá 127.921 ha
de erro fora da amostra para aquele modelo. A legenda passou a nomear a
especificação acoplada da Tabela S18a e a remeter à Tabela S35.

## 3. Números que divergiam do arquivo exportado

Todos na Tabela S42 e na S42b, e um deles também no corpo. Fonte:
`Outputs/Manuscrito/PartVI/TableS42_Sampling_Effort.csv` e
`TableS42b_Sampling_Effort_Robustness.csv`.

| Onde | Estava | Passou a ser | Valor de origem |
|---|---|---|---|
| S42, intercepto com esforço | 10,549 | **10,545** | 10,5449 |
| S42, IC do intercepto | [10,415, 10,684] | **[10,410, 10,680]** | 10,4098 e 10,68 |
| S42 e Seção 3.4, IC da forma | 0,665 | **0,666** | 0,666 |
| S42 e Seção 3.4, IC da proximidade | 0,619 | **0,620** | 0,6197 |
| S42, isolamento com esforço | 0,106 | **0,107** | 0,1067 |
| S42, IC do isolamento | [−0,045, 0,256] | **[−0,044, 0,257]** | −0,0439 e 0,2573 |
| S42, p do isolamento | 0,169 | **0,170** | 0,1698 |
| S42b, r de PLAND com log(OHE) | 0,130 | **0,129** | 0,1288 |
| S42b, jackknife da proximidade | 0,364 a 0,448 | **0,363 a 0,449** | 0,3632 e 0,4489 |
| S42b, jackknife da forma | 0,460 a 0,523, mediana 0,491 | **0,460 a 0,524, mediana 0,492** | 0,5238 e 0,4916 |
| S42b, proximidade sem unidades de zero registro | 0,397 [0,191, 0,603] | **[0,191, 0,604]** | 0,1906 e 0,6035 |
| S42b, forma sem unidades de zero registro | 0,456 [0,275, 0,636] | **[0,276, 0,637]** | 0,2755 e 0,6368 |
| S42b, proximidade sem as quatro mais amostradas | 0,403 [0,204, 0,601] | **[0,204, 0,602]** | 0,6020 |
| S42b, forma sem as quatro mais amostradas | 0,484 [0,303, 0,666] | **0,485 [0,303, 0,667]** | 0,4847 e 0,6669 |

O padrão é de truncamento no lugar de arredondamento, o que sugere que o corpo
dessas duas tabelas foi transcrito de uma execução anterior e só a legenda foi
revista na sessão de 01/09. As demais tabelas do suplementar não têm esse
padrão.

**Uma limitação que fica registrada.** O script 08 exporta com
`round(x, 4)`. Dois coeficientes caem exatamente no meio na terceira casa,
0,4125 para a proximidade e 0,4895 para a forma, e o arquivo de saída não tem
dígitos suficientes para decidir se a terceira casa é 2 ou 3, e 9 ou 0. Os
documentos trazem 0,412 e 0,489, que é o que estava, e as razões de extensão
declaradas, 1,51 e 1,63, são compatíveis com qualquer das duas leituras
(`OHE_ratio_per_SD` traz 1,511 e 1,632). Se essa casa vier a importar, basta
elevar a precisão do `round` nas linhas 231 a 234 do script 08 e reexportar.

## 4. Referências

Duas entradas constavam da lista do corpo sem serem citadas em lugar nenhum e
foram removidas: Legendre & Legendre (2012) e Solórzano et al. (2021b). Com a
saída da segunda, a citação e a entrada de Solórzano et al. (2021a) perderam o
sufixo, que existia só para desambiguar as duas.

Duas obras citadas na legenda da Tabela S18 não tinham entrada em nenhuma das
duas listas e foram acrescentadas à lista do suplementar, com o asterisco que
marca as obras citadas apenas ali: Aiken & West (1991) e Schielzeth (2010).
**Confira esses dois registros bibliográficos antes de submeter**; eles foram
escritos nesta sessão e não vieram de um arquivo do projeto.

## 5. Numeração e remissões

A Seção 2.5.2 não existia: o Método ia de 2.5.1 para 2.5.3. As três subseções
foram renumeradas (2.5.3 → 2.5.2, 2.5.4 → 2.5.3, 2.5.5 → 2.5.4) e as quatro
remissões cruzadas afetadas foram atualizadas, uma no corpo e três no
suplementar.

A legenda da Figura S10 remetia ao nulo de permutação da "Table S28a"; o nulo
está na S28b. Corrigido.

A introdução do suplementar citava "Tables S26 and S28", e não existe uma Tabela
S28 sem sufixo. Passou a "Tables S26 and S28a to S28f".

## 6. Precisão de linguagem

- A Seção 2.5.3 (antiga 2.5.4) chamava de "the untransformed response" o que é a
  raiz quadrada da área de classe absoluta, como o próprio
  `TableII_5_Response_Size_Dependence.csv` rotula. Reescrito, e as duas orações
  foram fundidas numa só.
- A Seção 3.2 e a legenda da Figura 4 chamavam a cobertura florestal de "measure
  of forest pattern". São seis medidas de configuração mais a cobertura, que é
  composição; e a cobertura é ajustada a uma ordenação construída sobre
  composição, o que explica o R² de 0,941. O texto passou a separar as duas
  coisas.
- A legenda da Figura 7 não explicava os asteriscos que aparecem no painel.
  Passou a trazer a convenção do script 03: \*\*\* para p < 0,001, \*\* para
  p < 0,01 e \* para p < 0,05.
- A Seção 4.1 repetia palavra por palavra a oração da Seção 2.5.1 sobre o fator
  de inflação da variância. A ocorrência da 4.1 foi encurtada.
- Legenda da Tabela S32: "does not exceed 3 × 10⁻⁹ in any column" vale para as
  16 colunas mostradas; entre as 35 comparadas, duas passam desse limite
  (LPI_Water 4,4 × 10⁻⁹ e LPI_Herbaceous 4,2 × 10⁻⁹). Passou a "in any column
  shown", e o expoente virou sobrescrito de verdade, em vez de `10^-9`.
- Legenda da Tabela S20 dizia "at least 0.998" onde o corpo diz "at or above
  0.999". Harmonizado com o corpo.
- Cabeçalho "2.5 Statistical Analysis " tinha espaço final e caixa alta fora do
  padrão das demais seções. Passou a "2.5 Statistical analysis".
- Texto de abertura da Tabela S1: `A two- stage`, `0. 85` (duas vezes),
  `large- sized`, `Table S 1` e `Non- vegetated` corrigidos.

## 7. Figuras

As figuras de origem estão todas a 600 dpi, com 4110 px na largura de 174 mm.
O suplementar embute os TIFFs originais nessa resolução. **O corpo não**: as oito
figuras estão embutidas a cerca de 220 dpi, e as Figuras 1 e 8, que vêm do
ArcGIS, a cerca de 150 dpi. Reinserir os TIFFs de `Outputs/Manuscrito/` no Word
resolve as seis figuras geradas em R; as Figuras 1 e 8 precisam ser reexportadas
do ArcGIS, item que já estava pendente.

A Figura 7 está recortada dentro do Word (5,7% à esquerda e 8,1% à direita) e a
Figura S1, verticalmente (17,9% acima e 11,3% abaixo). Verifiquei os dois
recortes contra o arquivo de origem: nenhum corta conteúdo, apenas margem.

Conferi o conteúdo das Figuras 6, 7 e S5 contra as tabelas. Os valores da Figura
7 batem célula a célula com a Tabela S24. As três estimativas da Figura 6 batem
com a Tabela 3. O painel da Figura S5 traz R² = 0,662, que é o valor da escala de
hectares, como agora diz a legenda.

Ponto menor: os eixos das figuras usam "Standardised", grafia britânica, e o
texto usa "standardized". Resolver exige alterar os scripts e reexportar.

## 8. Organização das pastas

Removido, com o repositório e a raiz do projeto voltando a ter uma cópia só de
cada coisa:

- `Repositorio_GitHub/Outputs/Manuscrito/`, 22 arquivos não rastreados gravados
  por uma execução do script 01 feita com o diretório de trabalho apontando para
  a pasta do repositório, às 12:57 de hoje. Os arquivos eram byte a byte
  idênticos aos da raiz. `Outputs/.gitkeep` permanece.
- `_Arquivo/Outputs_gravados_no_repo_29.08/` e `_Arquivo/Outputs_25.08_substituidos/`,
  saídas superadas pela execução de hoje.
- `_Arquivo/Conferencia_figuras_29.08/` e `_Arquivo/Docs_superados/`.
- Os três `.docx` superados de 29/08 e a versão anterior à correção da S42.
- `.Rhistory` da raiz do projeto.

`_Arquivo` passou de 18 MB para 6,8 MB e agora contém apenas as versões dos dois
documentos como estavam antes desta revisão, em
`_Arquivo/Manuscritos_superados/`.

Uma ressalva: `Repositorio_GitHub/Dados/Processados/Data_Raw_WithCoords.csv` foi
removido por engano no início da limpeza e restaurado com `git checkout`. Ele é
rastreado pelo Git e é a tabela derivada que a declaração de disponibilidade de
dados promete. O conteúdo restaurado é idêntico ao da raiz, diferindo apenas no
fim de linha.

## 9. O que continua pendente

**Bloqueia a submissão**

- `[REPOSITORY URL]`, `[DOI]` e `[VERSION]` na declaração de disponibilidade de
  dados, e os três campos correspondentes do `CITATION.cff`.
- Duas afiliações `[INSERT]` no cabeçalho e o telefone do autor correspondente.
- Publicar o repositório e depositar os 67 rasters no arquivo com DOI.

**Qualidade gráfica**

- Reinserir as seis figuras geradas em R no corpo, a 600 dpi.
- Reexportar as Figuras 1 e 8 do ArcGIS a 600 dpi.

**Reprodutibilidade**

- Rodar o script 01 a partir da raiz do projeto até o fim. A execução de hoje
  parou no passo [8/9]: `Dados/Processados/Data_Raw_WithCoords.csv` não foi
  reescrito, segue com data de 25/08, e falta `session_info_part01.txt` em
  `Outputs/Manuscrito/PartI`. Os 20 arquivos de saída da Parte I são byte a byte
  idênticos aos da execução completa das 12:57, e os dois
  `Data_Raw_WithCoords.csv` também, de modo que nenhum resultado depende disso;
  o que falta é o registro.

**Formatação, para decidir no Word**

- As legendas das Figuras 1 e 6 usam rótulo em negrito e texto em itálico; as
  das Figuras 2, 3, 4, 5, 7 e 8 estão inteiras em negrito e itálico. As Tabelas 1
  e 2 estão inteiras em negrito e a Tabela 3 só no rótulo.
- No suplementar, as Tabelas S36 e S37 aparecem depois da S42b e das figuras,
  porque estão na Seção 5. A numeração fica fora de ordem, ainda que a
  introdução avise onde elas estão.
- A lista de referências do suplementar não está inteiramente em ordem
  alfabética: Oliveira U e Tredennick estão entre Burnham e Dormann.

---

## Prompt desta sessão

```
Executei todos os scripts. Faça uma última conferência dos resultados
exportados, dos números e dos resultados relatados, das figuras e das legendas no
manuscrito e no material suplementar. Faça também uma revisão minuciosa,
cuidadosa e precisa do manuscrito e do material suplementar, comparando números e
tabelas. Se houver algo Redundante ou desnecessário, remova (seja nos scripts,
nos manuscritos e em outros documentos etc.). Me consulte sobre cada decisão
importante. Seja cuidadoso e detalhista. Não invente números e resultados. Evite
escrita robótica e padrões mecânicos e genéricos; seja claro, conciso e coerente.
A escrita deve ter uma lógica linear e prender o leitor, sem parágrafos nem
informações descontextualizadas. Tudo deve ser fundamentado na literatura
científica. Economize tokens sempre que possível.
```
