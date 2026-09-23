# Registro de correções — 15 de setembro de 2026

Sessão sobre a versão `Manuscript_Revisado_15.09.docx`, que chegou com três
comentários da revisão. Nenhum script foi executado; nenhum número foi
recalculado. Todo valor citado aqui foi lido dos CSVs de `Outputs/Manuscrito/`.

## 1. Comentário 1: o parágrafo do modelo a priori virou tabela

O parágrafo da Seção 3.4 que listava seis coeficientes em prosa passou a ser a
**Tabela 3**, com β, IC de 95%, p e a soma dos pesos de Akaike. As fontes são
`Table_S_PartIII_ConfirmatoryModel.csv` e
`Table_S_PartIII_VariableImportance.csv`. Incluí a linha do intercepto, que
estava na Tabela S18a mas não no parágrafo, para as duas tabelas de coeficientes
do corpo ficarem simétricas.

Sobraram do parágrafo três frases: o ponteiro para a tabela, os dois modelos
dentro de dois pontos de AICc que diferem só no isolamento, e a frase sobre o
produto negativo, que conduz ao parágrafo seguinte.

A tabela do modelo reportado passou a ser a **Tabela 4**. Renumerei as cinco
remissões do corpo e as catorze do suplementar, incluindo o cabeçalho de coluna
da Tabela S42. As Tabelas S18a e S18b continuam no suplementar, com o erro-padrão
e as médias de modelo, e a legenda da nova Tabela 3 remete a elas.

Cópia em Excel: `Manuscrito/Manuscritos_Word/Table3_Full_APriori_Model.xlsx`.

## 2. A passagem do teste de permutação

A revisão de 15.09 reescreveu a Seção 3.4 e trocou o que o teste mede. O texto
passou a dizer que ele avalia "which coefficients are statistically significant",
que os coeficientes "are not statistically significant", e que os percentis 2,8,
6,9 e 10,4 são "what would be expected if no real effect existed". As três
afirmações contradizem o argumento do artigo: a Tabela S28c mostra justamente que
o critério convencional não arbitra aqui, e sob a nula os percentis são uniformes,
de modo que 2,8 não é mais esperado que qualquer outro.

Reescrevi as duas passagens. O teste agora é descrito pelo que faz — embaralha os
numeradores das duas métricas normalizadas e mantém cada denominador no valor
observado — e o resultado pelo que ele é: os coeficientes caem dentro do intervalo
central de 95% da distribuição que a construção da razão gera sozinha. Os seis
números citados não mudaram.

Também corrigi a frase truncada "since it is two-sided below the 5% significance
level" e, na Seção 4.1, "prevents this", que se referia à separação e não ao
controle.

## 3. Comentário 3: a Seção 4.6

O pedido era apagar a seção inteira. Ela continha dois resultados citados em
outras partes do artigo, e não apenas ressalvas. A solução foi separar as duas
coisas.

Subiram para o fim da Seção 4.1, onde são argumento e não defeito: o raio de busca
do índice de proximidade, 28.230 m, que excede o diâmetro circular equivalente de
54 das 67 unidades; o desbaste espacial, que reduz o desenho de 67 para 56 unidades
e, em 200 réplicas, deixa as duas associações em pé e ligeiramente mais fortes,
com medianas de 0,81 e 0,66 (Tabela S41); e a dependência do esforço amostral,
cujas estimativas condicionais são um piso e não uma correção (Tabela S42).

A seção passou a se chamar **4.6 Scope and further work**, com três parágrafos: o
alcance dos resultados, o caráter estático da descrição, e o que vem depois. Sem
essa mudança, a Tabela S41 ficaria sem citação no corpo e o raio de busca
desapareceria do texto.

## 4. Comentário 2: a Seção 4.5

O parágrafo longo sobre os limites do Índice de Valor de Conservação foi comprimido
a duas frases que remetem à Seção 5 do suplementar, onde os limites já estão por
extenso. As quatro recomendações perderam o "First, Second, Third, Fourth" e o
imperativo. O parágrafo que rotula duas leituras como hipóteses, e não como
recomendações, ficou.

Também corrigi a abertura da seção, que na versão recebida dizia "based on the two
associations that are not the coupling diagnostic and spatial control".

## 5. Uma frase removida da Seção 4.3

A versão de 15.09 acrescentara "This paper suggests that sparing is the optimal
conservation strategy for the state". Nenhuma análise do artigo compara as duas
estratégias; o que há é uma associação transversal, o envelope 1,49 vez maior por
desvio padrão de consolidação da matriz (Tabela S40), que a própria Seção 4.5
classifica como hipótese a testar. A frase também colidia com a recomendação de
agrofloresta logo em seguida. Restou a oração que já estava lá: a separação estrita
é raramente viável no estado.

## 6. Um ponto a confirmar

A Seção 2.2 define porte médio como 1,0 a 7 kg, citando Chiarello (2000), e a
Seção 4.6 acompanha. Até 02 de setembro o limiar era 1,5 kg. Confirme que
Chiarello (2000) sustenta 1,0 kg; é a única fonte dada para o limiar.

## 7. Estado dos arquivos

Versões anteriores em `_Arquivo/Manuscritos_superados/`, com os sufixos
`_antes_tabela` e `_antes_revisao_discussao`. Os três comentários da revisão
continuam no arquivo, com as marcas de âncora intactas, para serem resolvidos no
Word. As oito imagens seguem byte a byte idênticas às da versão recebida.
