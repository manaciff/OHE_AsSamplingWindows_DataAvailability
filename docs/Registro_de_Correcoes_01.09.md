# Registro de correções — 01 de setembro de 2026

Sessão dedicada a três coisas: verificar se a ambiguidade do teste de permutação
apontada na auditoria de 22 de agosto ainda existia, rastrear os números do corpo
e do suplementar até a célula que os produz, e deixar o repositório em condição de
ser publicado. Todo valor citado aqui foi lido do arquivo de saída indicado.

---

## 1. O teste de permutação: a ambiguidade não existe mais

A auditoria de 22 de agosto registrou que os valores da nula eram chamados de
percentis quando eram valores de p empíricos, e que duas regras de decisão
conviviam sem que o texto dissesse qual valia. Conferi os dois documentos.

O corpo, na Seção 3.4, diz:

> The observed coefficients for forest cover, patch density, and the product
> returned two-sided empirical p-values of 0.029, 0.075, and 0.155 against the
> null, and they sit at the 2.8th, 6.9th, and 10.4th percentiles of their
> respective null distributions. All three fall within the central 95% interval
> of their null distributions, which is the criterion we apply, because a
> comparison of magnitudes alone cannot summarize a null that is neither
> symmetric nor centered on zero.

O suplementar, na legenda da Tabela S28b, repete a regra e trata o caso limítrofe:

> Forest cover is the borderline case: it sits at the 2.8th percentile of its null
> distribution, so the two-sided empirical p is 0.0294, below 0.05, while the
> coefficient itself remains inside the interval. The two readings are reported
> side by side rather than reconciled.

As duas quantidades estão nomeadas corretamente, a regra está declarada, a
justificativa está dada e o caso limítrofe está reconhecido. **Nada a ajustar.**

Os seis valores conferem com `PartV/TableS_PartV_NullSimulation.csv`:

| Termo | Percentil | p bilateral |
|---|---|---|
| PLAND_Forest | 2.81 | 0.0294 |
| PD_Forest | 6.91 | 0.0748 |
| PLAND × PROX | 10.38 | 0.1545 |

### Uma correção na nota que circulou

O texto do orientador que motivou esta verificação lista os percentis como
"2,9%, 7,5% e 15,5%". Os percentis são **2,8%, 6,9% e 10,4%**. O 15,5% é o valor
de p da interação, 0,1545, e não o percentil dela. Vale corrigir a nota antes que
ela seja reutilizada.

---

## 2. As duas correções aplicadas

Divergência real entre os dois documentos, na legenda da Tabela S42 do suplementar.
Fonte: `PartVI/TableS42_Sampling_Effort.csv`.

| Onde | Estava | Passou a ser | Valor de origem |
|---|---|---|---|
| Legenda S42 | deviance de 36.1% | **36.0%** | 0.3604, modelo M1 |
| Legenda S42 | sobe para 57.0% | **56.9%** | 0.5692, modelo M2 |

O corpo já dizia 36.0%, de modo que os dois documentos se contradiziam. Os fatores
de inflação da variância citados na mesma frase, 1.93 e 2.32, estão certos:
o CSV traz 1.934 e 2.316.

O arquivo anterior está em
`_Arquivo/Manuscritos_superados/Supplementary_Material_REVISED_30.08_antes_correcao_S42.docx`.

---

## 3. O que foi conferido e está correto

Rastreado até o arquivo de origem, não apenas localizado no conjunto de valores
exportados.

| Afirmação | Valor | Fonte |
|---|---|---|
| Mediana de cobertura florestal | 54.9% | 54.8547, `Table_S_PartIII_DescriptiveStatistics.csv` |
| Unidades abaixo de 30% | 18 | contagem direta em `Data_Raw_WithCoords.csv` |
| Razão por desvio padrão, proximidade | 2.10 | 2.1, `TableS_PartV_Uncoupled_Model.csv` |
| Razão por desvio padrão, forma | 1.83 | 1.826, mesmo arquivo |
| Razões controladas por esforço | 1.51 e 1.63 | 1.511 e 1.632, `TableS42_Sampling_Effort.csv` |
| Matriz agropecuária consolidada | 1.49 | 1.492, `TableS_PartV_Matrix_Uncoupled.csv` |
| Densidade de manchas em *Brachyteles* | 0.63 | corrigido de 0.64 em revisão anterior |
| Deviance do modelo acoplado e do reportado | 55.1% e 36.0% | atribuição explícita no corpo |
| Erro de validação cruzada | raiz do erro quadrático médio | nomeado corretamente |
| Coeficiente de agrupamento | exclui zero em cinco das seis especificações | alinhado com a Tabela S28a |
| Amplitude da resposta | 1.906,73 a 386.509,44 ha | precisão total, sem arredondamento |

O último item merece nota. A auditoria de 22 de agosto tinha registrado que a
variável resposta estava gravada como inteiro, com 28 dos 67 valores múltiplos de
100. Isso foi resolvido: `arcpy_AOH_area_full_precision.py` reexportou a área com
precisão total, e num artigo que discute uma identidade algébrica entre resposta e
preditores essa precisão é parte do argumento.

---

## 4. Itens da auditoria de 22 de agosto já resolvidos

Verificados no código e no texto, não presumidos.

- **M1, a ordenação canônica.** O script 02, linha 924, aplica
  `vegan::decostand(LULC_CA, method = "hellinger")`. O corpo descreve a RDA sobre
  a área de classe transformada por Hellinger, e a Tabela S39 documenta o efeito
  da transformação.
- **M2, as frações da partição.** O corpo diz agora que 9,9% foi compartilhado com
  o espaço e que a fração espacial pura ficou em 1,2%, que é a ordem coerente com
  o teste de permutação daquela fração.
- **M3, a métrica acoplada na partição.** A Tabela S38 reporta a partição sem
  `PD_Forest`.
- **M5, os testes por termo.** O script 02, linha 993, usa `by = "margin"`, e o
  corpo declara que reporta testes marginais.
- **Entradas alternativas do script 03.** O script exige `Data_Raw_WithCoords.csv`
  e interrompe a execução com uma mensagem explicativa se ele faltar.
- **Atribuição de Rubia Morini.** Consta apenas dos agradecimentos; a lista de
  contribuições traz só os cinco autores.

---

## 5. O repositório

Estado inicial: pasta chamada `Repositorio_GitHub` sem `git init`, portanto sem
histórico e sem como publicar ou gerar DOI.

O que mudou:

1. **`git init` e commit inicial.** 138 arquivos, 12 MB. O `.gitignore` já estava
   correto e manteve fora os 67 rasters de 2,6 GB, distribuídos por arquivo com
   DOI, e os arquivos de sessão.
2. **`R/run_demo3.R` removido.** Chamava `source('/home/claude/fix/R/06_connectivity_demo.R')`,
   caminho de uma sessão antiga. Verifiquei: os scripts 06 e 07 e as duas saídas
   da demonstração de conectividade não existem em nenhum lugar do projeto. Era o
   último fragmento de um conjunto já incompleto e não produzia nada citado no
   manuscrito.
3. **`.RData` de 10 MB e `.Rhistory` removidos do disco.** Estavam no `.gitignore`,
   mas o `.RData` carrega objetos ao abrir o R na pasta e pode injetar resultados
   de execuções anteriores numa sessão sem aviso.
4. **Script 08 incluído no `run_all.R`.** Ele escreve as Tabelas S42 e S42b, que o
   manuscrito cita. Quem clonasse o repositório e rodasse `run_all.R` não as
   reproduzia. Roda depois do 05, porque lê o modelo da Tabela 3.
5. **`docs/SCRIPTS.md`.** Ganhou a seção do script 08. A seção sobre a
   demonstração de conectividade foi reescrita: o argumento sobre PC e dPC ficou,
   porque conversa diretamente com a Seção 2.5.1 do manuscrito, e o texto agora
   declara que os scripts e as figuras não estão arquivados aqui.
6. **`docs/TABLE_MAP.md`.** Ganhou as linhas das Tabelas S42 e S42b, que faltavam
   e contradiziam a promessa do próprio arquivo de rastrear todo número reportado.
   A seção de conectividade foi reescrita pelo mesmo motivo do item anterior.
7. **`README.md`.** A linha que descrevia o `run_all.R` dizia "runs the five
   scripts in order" sem explicar por que os scripts 00 e 01 estão comentados.

---

## 6. O que continua pendente

Ordenado por quem bloqueia o quê.

**Bloqueia a submissão**

- Três marcadores na declaração de disponibilidade de dados do corpo,
  `[REPOSITORY URL]`, `[DOI]` e `[VERSION]`, e os três correspondentes no
  `CITATION.cff`: `date-released`, `repository-code` e `doi`. Resolvem-se juntos,
  quando o repositório remoto existir e o depósito no Zenodo gerar o DOI.
- Duas afiliações marcadas com `[INSERT]` no cabeçalho.
- Publicar o repositório e depositar os 67 rasters no arquivo com DOI.

**Bloqueia a reprodução completa**

- Rodar o script 01 até o fim. Em `Outputs/Manuscrito/PartI` há 21 arquivos, 12 de
  29 de agosto e 9 ainda de 25 de agosto, o que mostra que a última execução parou
  no meio. Os valores conferem onde há sobreposição, mas a execução precisa ser
  completa antes da defesa.

**Qualidade gráfica**

- Reexportar as Figuras 1 e 8 do ArcGIS a 600 dpi. As demais já estão nessa
  resolução.

---

## 7. O prompt desta sessão

```
verifique e ajuste esse ponto se a ambiguidade de fato existir no manuscrito:
"4. Teste de Permutação Nula como Filtro de Artefato
Status: plenamente implementado, Parte V
Seu teste em 9.999 permutações prova que:
 - Proximidade florestal: 2,9% de p empírico (dentro do intervalo central de 95%
   da nula, não distinguível por esse critério)
 - Complexidade de forma: 7,5% (idem)
 - Produto: 15,5% (idem)
Mas há uma ambiguidade que a auditoria flagou (M4): você reporta esses como
'percentis' quando são valores de p empíricos. A regra de decisão adotada é o
intervalo central de 95%, não o valor de p. Essas duas regras discordam para
PLAND, então você precisa declará-las explicitamente."

Organize novamente a pasta e o repositório do GitHub de forma precisa e detalhada.
Faça também uma revisão minuciosa, cuidadosa e precisa no manuscrito e no material
suplementar, comparando números e tabelas. Se houver algo redundante ou
desnecessário, remova (seja nos scripts, nos manuscritos e em outros documentos
etc.). Me consulte sobre cada decisão importante. Seja cuidadoso e detalhista.
Não invente números e resultados. Evite escrita robótica e padrões mecânicos e
genéricos; seja claro, conciso e coerente. A escrita deve ter uma lógica linear e
prender o leitor, sem parágrafos ou informações descontextualizadas. Tudo deve ser
fundamentado na literatura científica. Economize tokens sempre que possível. Ao
final faça um arquivo md com o resumo dessa sessão e com esse prompt salvo/escrito
no md.
```

---

## Nota sobre redundância

Rodei uma comparação par a par das frases do corpo e do suplementar. O que
apareceu foram legendas paralelas por construção, as PERMANOVA par a par por
gênero, espécie, dieta e locomoção, que descrevem tabelas diferentes e cuja
uniformidade ajuda a leitura, e entradas da lista de referências de mesmos
autores. Dois pares no corpo repetem informação entre a Seção 2.5.5 e a Seção 3.4,
o tamanho amostral por gênero e o anúncio do teste de permutação, mas cada
ocorrência faz um trabalho diferente: uma justifica uma decisão de desenho e a
outra qualifica a leitura do resultado. Não removi nenhuma. Não encontrei
parágrafo dispensável nos dois documentos.
