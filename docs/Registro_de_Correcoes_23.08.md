# Registro de correções, 23 de agosto de 2026

Sessão posterior à execução do `arcpy_AOH_area_full_precision.py` no ArcGIS e à
reexecução dos seis scripts. Todos os números citados abaixo vêm dos CSV
exportados pela cadeia; nenhum foi digitado à mão.

---

## 1. Defeitos encontrados nesta sessão

### 1.1 O script 04 lia um arquivo diferente dos scripts 03 e 05

`04_geometric_coupling_diagnostic.R` lia `Data_Raw_FINAL.csv`, no qual a resposta
está gravada como inteiro, enquanto os scripts 03 e 05 leem
`Data_Raw_WithCoords.csv`, que carrega a resposta com precisão total. As duas
tabelas do suplementar que verificam a mesma identidade algébrica traziam, por
isso, valores diferentes: a Tabela S26 reportava discrepância máxima de 0,0290 e a
Tabela S28e reportava 0,0277 para a mesma quantidade, nas mesmas 67 unidades.

**Consequência.** Um revisor que confrontasse as duas tabelas encontraria uma
inconsistência interna sem explicação no texto.

**Correção.** O script 04 passa a ler `Data_Raw_WithCoords.csv`, com a mesma
sequência de nomes alternativos de coluna que o script 03 usa
(`AOH_unit_area_ha`, `UA_Area_ha`, `EHA_ha`). O guard de média de `log(OHE)` foi
atualizado para 10,3735 com tolerância de 0,01, que aceita tanto a resposta
arredondada quanto a de precisão total e rejeita a versão com o erro de junção.
A Tabela S26 foi refeita com esse arquivo. A identidade agora fecha em 0,0277 nas
duas tabelas.

### 1.2 As legendas das Figuras S5 e S6 descreviam outro modelo

A Seção 8 do script 03 escolhe a especificação retida entre os modelos dentro de
duas unidades de AICc, primeiro pelo menor número de parâmetros e depois pelo
menor erro de validação cruzada. Nestes dados isso seleciona o modelo Gamma com
ligação logarítmica, de modo que `global_model` é um `glm` e todas as figuras
seguintes são desenhadas a partir dele. As duas legendas descreviam um modelo
misto, previsão em nível populacional e um quarto painel de resíduos por gênero,
nada disso presente nas figuras.

A conferência é direta: o painel da Figura S5 traz R² = 0,662, que é a correlação
ao quadrado entre observado e predito na escala de hectares para o modelo Gamma,
e a Figura S6 mostra resíduos de Pearson, resíduos de deviance e um painel de
resíduos contra alavancagem, que um modelo misto não produz aqui.

**Correção.** As duas legendas foram reescritas. A da Figura S5 agora reporta
0,662 na escala de hectares e 0,557 na escala logarítmica, e a maior unidade
amostral com o valor correto de 386.509 ha. Os três valores que a legenda antiga
citava, 0,43, 0,54 e 0,13, não correspondem a nenhum objeto da cadeia atual e
foram retirados, junto com a remissão a colunas de R² marginal e condicional que
a Tabela S16 não tem.

**Prevenção.** O script 03 passa a exportar `Table_S_PartIII_ObsPred_R2.csv` com
as duas correlações ao quadrado, o rótulo do modelo usado e a maior unidade
amostral, de modo que a legenda deixa de depender de leitura da imagem.

### 1.3 A Tabela S19 tinha duas colunas duplicadas

O cabeçalho e as duas linhas de dados traziam ΔAICc e w duas vezes, em nove
colunas. As colunas repetidas foram removidas e as larguras redistribuídas, de
modo que os valores numéricos deixaram de quebrar em duas linhas.

### 1.4 Um limite superior arredondado para baixo

A legenda da Tabela S34 afirmava que nenhum coeficiente se moveu mais que 0,16. A
maior mudança absoluta é 0,1633, de modo que o limite era falso. Passou a 0,17,
que é o valor arredondado para cima e coincide com o que o manuscrito já
reportava nas Seções 3.4 e 4.1.

### 1.5 Três pacotes listados e não citados

DHARMa, car e MuMIn constavam da lista de referências do suplementar sem chamada
no texto. Foram citados onde cada um é usado: DHARMa na legenda da Tabela S16,
car na da Tabela S18 e MuMIn na da Tabela S19.

### 1.6 Remissões imprecisas

A legenda da Tabela S23 remetia à "Tabela S28", que não existe como tabela única,
e passou a remeter às Tabelas S28a a S28e. A legenda da Tabela S38 remetia às
frações da "Tabela S27", que é a tabela de tamanhos de amostra, e passou a
remeter à Seção 3.3 do manuscrito, onde as frações são de fato reportadas.

### 1.7 O `run_all.R` não encontrava os scripts

Com `ANALISES_TCC.Rproj` aberto, `here()` resolve para `D:/Duda_Nacif_TCC`, onde
não existe a pasta `R/`, e a linha `source(here("R", s))` falhava na primeira
iteração. O `run_all.R` passa a localizar a cadeia a partir da própria posição do
arquivo, com uma lista curta de pastas alternativas quando isso não é possível, e
`here()` continua resolvendo `Dados/` e `Outputs/` dentro da pasta de trabalho.
As duas pastas são verificadas antes do primeiro `source`, com mensagem
explicativa se faltarem.

### 1.8 Nome da resposta

As tabelas exportadas pelos scripts 04 e 05 rotulavam a resposta como `log(AOH)`,
enquanto o manuscrito a chama de OHE, envelope de habitat ocupado. Os rótulos
exportados e os rótulos dos eixos das figuras do script 03 passaram a usar OHE.
No suplementar, quatorze rótulos de célula foram renomeados.

---

## 2. O que mudou nos resultados com a resposta de precisão total

A resposta deixou de estar arredondada para a dezena, a centena ou o milhar mais
próximos. O efeito é pequeno e nenhuma conclusão muda.

| Quantidade | Antes | Agora |
|---|---|---|
| Média de OHE (ha) | 47.836,33 | 47.639,02 |
| Maior unidade (ha) | 387.619 | 386.509 |
| Menor unidade (ha) | 1.910 | 1.906,73 |
| Identidade algébrica, discrepância máxima | 0,0290 | 0,0277 |
| Razão de extensão por desvio padrão de proximidade | 2,10 | 2,10 |
| Razão de extensão por desvio padrão de complexidade de forma | 1,82 | 1,83 |
| Coeficiente de proximidade, modelo reportado | 0,7413 | 0,7420 |
| Fração puramente espacial, conjunto completo | 1,2% | 1,3% |
| Fração puramente espacial, sem PD | 6,7% | 6,8% |
| Permutações que alcançam a vantagem de AICc | 38,8% | 38,1% |
| RMSE de validação cruzada, modelo Gamma | 40.102 | 39.977 |

A razão da complexidade de forma volta a 1,83 porque exp(0,6023) = 1,8263. Com a
resposta arredondada o coeficiente era 0,6015 e a razão saía 1,82.

---

## 3. O que foi reescrito no material suplementar

Vinte e sete tabelas foram reescritas célula a célula a partir dos CSV
correspondentes, por um script que lê o rótulo de cada linha do documento,
procura o registro correspondente no CSV e formata o valor por uma regra
declarada no código: S15, S16, S17, S18a, S18b, S19, S20, S21, S23, S24, S25,
S26, S28a, S28b, S28c, S28d, S28e, S28f, S29, S30, S31, S34, S35, S38, S40, S41.
A Tabela S30 não teve nenhuma célula alterada, o que confirma que a regra de
arredondamento adotada reproduz a que estava em uso.

Não mudaram: S1 a S14 e S27, S32, S33, S39 e os Apêndices S1 a S3. A Parte I e a
ordenação da Parte II não dependem da área da unidade; a Tabela S39 foi
reconferida linha a linha contra
`TableII_4b_RDA_Transformation_Comparison.csv` e está idêntica.

Além das tabelas, foram atualizadas as legendas de S17, S20, S21, S28b, S28c,
S34, S38 e S40 e as legendas das Figuras S5 e S6.

A legenda da Tabela S28b passou a declarar o caso limítrofe: a cobertura florestal
está no percentil 2,8 da sua distribuição nula, de modo que o p empírico
bilateral é 0,0294, abaixo de 0,05, enquanto o coeficiente permanece dentro do
intervalo central de 95%. As duas leituras são reportadas lado a lado, e a
conclusão que a tabela sustenta, que o sinal desse coeficiente não deve ser
interpretado, é a mesma nas duas.

---

## 4. Verificações executadas

1. Validação XSD dos dois documentos contra os originais: aprovada. O suplementar
   perdeu seis parágrafos, que são exatamente as seis células removidas das duas
   colunas duplicadas da Tabela S19.
2. Conferência numérica automática: cada número impresso nos dois documentos foi
   procurado no conjunto de todos os CSV exportados, em todos os arredondamentos
   de zero a seis casas. Sobraram apenas anos de citação, DOI, números de página,
   valores metodológicos declarados no texto, como 28.230 m e 42,4 m, e os
   números do índice de valor de conservação, que vêm do ArcGIS e não da cadeia
   em R.
3. Itálico: nenhuma ocorrência de nome de gênero ou de espécie fora de itálico,
   nos dois documentos.
4. Referências: no manuscrito, 107 entradas, 105 citadas. As duas restantes,
   Magioli et al. (2015) e Püttker et al. (2020), continuam marcadas em amarelo
   aguardando decisão. No suplementar, 18 entradas, todas citadas.
5. Reprodução independente do modelo Gamma sobre `Data_Raw_WithCoords.csv`: os
   sete coeficientes reproduzem `Table_S_PartIII_ConfirmatoryModel.csv` com
   quatro casas decimais, e o R² na escala de hectares reproduz o valor impresso
   dentro da Figura S5. Foi essa reprodução que permitiu corrigir as legendas sem
   inventar valores.
6. Os oito arquivos `.R` foram lidos por `parse()` sem erro.

---

## 5. Pendências

1. Regerar as Figuras S5 e S6 após a próxima execução do script 03. A S5 porque
   os rótulos dos eixos passaram de AOH para OHE, e a S6 porque o título do painel
   estava sobreposto ao texto que `plot.lm` escreve na margem externa, agora
   suprimido com `sub.caption = ""`.
2. Decidir sobre Magioli et al. (2015) e Püttker et al. (2020): citar no corpo ou
   retirar da lista.
3. Preencher os três marcadores de `CITATION.cff`: `date-released`,
   `repository-code` e `doi`.
4. Preencher as duas afiliações marcadas em amarelo na folha de rosto do
   manuscrito.
