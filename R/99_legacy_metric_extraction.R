# ==============================================================================
# 99  Legacy metric extraction with landscapemetrics (FRAGSTATS-like)
#
# Manuscript: Assessing the Landscape Composition and Configuration of
#             Threatened Medium- and Large-Sized Mammals using Habitat
#             Envelopes as Sampling Windows
# Author:     Maria Eduarda Nacif
#
# Superseded by script 00. Kept because it documents the provenance of the
# consolidated dataset and of the class aggregation. Do not run it as part of
# the current chain.
#
# The guard below enforces that. Until 24 August 2026 the header said the file
# must not be run while the file itself was fully executable and wrote three CSV
# into Dados/Processados, one of which, Data_Raw_GLM_NMDS_Final.csv, script 03
# names as a file that yields a different sample. Sourcing this file by accident
# overwrote inputs of the chain. To read it as documentation, open it. To run it
# deliberately, set RUN_LEGACY_99 <- TRUE first.
# ==============================================================================
if (!isTRUE(getOption("run.legacy.99", FALSE)) &&
    !isTRUE(mget("RUN_LEGACY_99", envir = globalenv(),
                 ifnotfound = list(FALSE))$RUN_LEGACY_99)) {
  stop("99_legacy_metric_extraction.R is superseded documentation and writes ",
       "over inputs of the current chain. Set RUN_LEGACY_99 <- TRUE before ",
       "sourcing it if that is really what you want.", call. = FALSE)
}
# ==============================================================================
#
# PROCEDÊNCIA DO MAPA DE COBERTURA (documentado em 29/07/2026)
#
# Fonte........: MapBiomas, Coleção 10.1 (série 1985-2024), resolução 30 m.
# Ano de refer.: 2024 (recorte único; ver Dados/Rasters/Lulc_2024.tif).
#                Os rasters por unidade em UA_RASTER/ derivam desse mesmo ano.
# Acurácia.....: Mata Atlântica, 91,5% no nível 1 e 86,1% no nível 2, validada
#                com 10.840 pontos independentes do LAPIG/UFG.
#                Fonte: Atlantic Forest Appendix, ATBD Collection 10, v1.
# UMM..........: aproximadamente 0,5 ha, imposta por filtro espacial de seis
#                pixels conectados, com exceção que preserva floresta ripária.
#                Esse limiar define o menor fragmento que pode entrar em NP e
#                PD, e por isso precisa constar dos Métodos.
#
# ATENÇÃO À CONSISTÊNCIA TEMPORAL: os registros de ocorrência que delimitam as
# AOH cobrem 1990-2020, enquanto a paisagem medida aqui é a de 2024. O manuscrito
# ainda não declara o ano de referência do mapa. É a informação que o parecer
# pede em "temporal consistency" e precisa entrar na seção 2.4.
#
# As classes 1..6 usadas na matriz de reclassificação abaixo já são o produto de
# uma agregação prévia das 22 subclasses terminais do MapBiomas (ver Tabela S1).
# Aqui a Agricultura (4) é ainda fundida na Pastagem (3), formando a categoria
# "Agropecuária" (Farming) usada no manuscrito.
# ==============================================================================

# 1. Carregar Pacotes
library(terra)            # Manipulação de rasters
library(landscapemetrics) # Cálculo das métricas (Substituto do FRAGSTATS)
library(tidyverse)        # Manipulação de dados (dplyr, tidyr)
library(here)             # Caminhos relativos
#citation("landscapemetrics")

# Caminhos (relativos ao .Rproj na raiz D:/Duda_Nacif_TCC)
data_path  <- here("Dados", "Processados")  # CSVs intermediários (input/output)
rasters_path <- here("Dados", "Rasters")    # Rasters de entrada
fragstats_path <- here("Dados", "FRAGSTATS_RESULT")

# 2. Definir a pasta dos rasters

pasta_rasters <- file.path(rasters_path, "UA_RASTER")
lista_rasters <- list.files(pasta_rasters, pattern = "\\.tif$", full.names = TRUE)

message("Encontrados ", length(lista_rasters), " rasters na pasta.")

# 3. Matriz de Reclassificação (Agrupando Agricultura [4] na Pastagem [3])
reclass_matrix <- matrix(c(
  1, 1,  # Floresta
  2, 2,  # Herbácea e Arbustiva
  3, 3,  # Pastagem 
  4, 3,  # Agricultura -> Vira classe 3 (Agropecuária)
  5, 5,  # Água
  6, 6   # Área Não Vegetada
), ncol = 2, byrow = TRUE)

# 4. Loop para processar cada raster
resultados_lista <- list()

for(ficheiro in lista_rasters) {
  
  nome_base <- basename(ficheiro)
  nome_sem_ext <- tools::file_path_sans_ext(nome_base)
  
  partes <- strsplit(nome_sem_ext, "_")[[1]]
  ud_id <- partes[1]
  especie <- partes[2]
  
  rst <- rast(ficheiro)
  rst_reclass <- classify(rst, reclass_matrix)
  
  # CÁLCULO DAS MÉTRICAS DE CLASSE (Sem o PROX)
  metricas <- calculate_lsm(rst_reclass, 
                            level = "class", 
                            metric = c("ca", "pland", "np", "pd", "area_mn", 
                                       "enn_mn", "frac_mn", "lpi", "ed"), 
                            directions = 8)
  
  metricas_limpas <- metricas %>%
    select(class, metric, value) %>%
    mutate(
      UA_ID = ud_id,
      Species_Code = especie,
      Class_Name = case_when(
        class == 1 ~ "Forest",
        class == 2 ~ "Herbaceous",
        class == 3 ~ "Agropecuaria", 
        class == 5 ~ "Water",
        class == 6 ~ "Non_Vegetated",
        TRUE ~ as.character(class)
      )
    )
  
  resultados_lista[[nome_sem_ext]] <- metricas_limpas
}

# 5. Juntar e pivotar
tabela_final_long <- bind_rows(resultados_lista)

tabela_final_wide <- tabela_final_long %>%
  pivot_wider(
    id_cols = c(UA_ID, Species_Code),
    names_from = c(metric, Class_Name),
    values_from = value,
    names_glue = "{toupper(metric)}_{Class_Name}" 
  )

write_csv(tabela_final_wide, file.path(data_path, "Tabela_Metricas_Landscapemetrics.csv"))
message("\n >> Sucesso! Tabela salva.")

# ==============================================================================
# SCRIPT 2: RESGATANDO O PROX DO FRAGSTATS E INTEGRANDO AOS DADOS
# ==============================================================================

library(tidyverse)
library(here)

message("--- Iniciando o resgate do PROX_MN ---")

# 1. Encontrar todos os arquivos de Patch na pasta FRAGSTATS_RESULT
pasta_fragstats <- fragstats_path
# Procura por arquivos que terminem com "Patch.csv" (ignora maiúsculas/minúsculas)
arquivos_patch <- list.files(pasta_fragstats, pattern = "Patch\\.csv$", ignore.case = TRUE, full.names = TRUE)

message("Encontrados ", length(arquivos_patch), " arquivos do FRAGSTATS.")

# 2. Ler e juntar todos os arquivos num só dataframe (À Prova de Falhas)
dados_fragstats_patch <- arquivos_patch %>%
  map_dfr(~ {
    # Ler o arquivo
    df <- read_csv(.x, col_types = cols(.default = "c"), name_repair = "unique_quiet")
    nome_arquivo <- basename(.x)
    
    # Checar se a coluna PROX existe. Se não existir, avisa e cria uma vazia!
    if (!"PROX" %in% names(df)) {
      message(" ⚠️ AVISO: O arquivo '", nome_arquivo, "' NÃO possui a coluna PROX. Preenchendo com NA.")
      df$PROX <- NA
    }
    
    # Pegamos apenas as colunas que interessam de forma segura
    df %>%
      select(any_of(c("SPECIES", "UA", "TYPE", "PROX"))) %>%
      # Forçamos o PROX a voltar a ser numérico
      mutate(PROX = as.numeric(PROX))
  })

# 3. Calcular o PROX_MN por Classe e Unidade Amostral
prox_compilado <- dados_fragstats_patch %>%
  # Remover espaços em branco que o FRAGSTATS costuma deixar (ex: " cls_1 ")
  mutate(TYPE = str_squish(TYPE)) %>%
  # Criar a nossa coluna de nomes amigáveis para bater com a tabela do R
  mutate(Class_Name = case_when(
    TYPE == "cls_1" ~ "Forest",
    TYPE == "cls_2" ~ "Herbaceous",
    TYPE %in% c("cls_3", "cls_4") ~ "Agropecuaria", # Agrupa pasto e agricultura
    TYPE == "cls_5" ~ "Water",
    TYPE == "cls_6" ~ "Non_Vegetated",
    TRUE ~ TYPE
  )) %>%
  # Agrupar por Espécie, UA e a nova Classe
  group_by(SPECIES, UA, Class_Name) %>%
  # Calcular a média do PROX (ignorando manchas que possam ter NA)
  summarise(PROX_MN = mean(PROX, na.rm = TRUE), .groups = "drop") %>%
  # Pivotar para o formato largo (Wide)
  pivot_wider(
    id_cols = c(UA, SPECIES),
    names_from = Class_Name,
    values_from = PROX_MN,
    names_prefix = "PROX_MN_"
  )

# 4. Carregar a nossa tabela gerada pelo landscapemetrics
tabela_landscapemetrics <- read_csv(file.path(data_path, "Tabela_Metricas_Landscapemetrics.csv"), show_col_types = FALSE)

# IMPORTANTE: No FRAGSTATS, a sua coluna de espécie pode estar inteira (ex: "M_tridactyla"), 
# enquanto no raster estava abreviada ("tridactyla1"). 
# Se os nomes na coluna SPECIES não baterem perfeitamente com a Species_Code, 
# precisaremos de um ajuste rápido aqui antes do Join. Mas vamos tentar juntar pelo UA_ID.

# Ajustando o nome da coluna UA para o Join
prox_compilado <- prox_compilado %>% rename(UA_ID = UA)

# 5. Juntar as duas tabelas! (Left Join)
# Vamos juntar usando o UA_ID (Unidade Amostral)
tabela_final_absoluta <- tabela_landscapemetrics %>%
  # Garante que o ID seja do mesmo tipo (texto)
  mutate(UA_ID = as.character(UA_ID)) %>%
  left_join(
    prox_compilado %>% mutate(UA_ID = as.character(UA_ID)), 
    # Juntamos pelo ID da Unidade Amostral e pelo Código da Espécie
    # Se os códigos de espécie não baterem, avise-me!
    by = c("UA_ID" = "UA_ID", "Species_Code" = "SPECIES") 
  )

# Salvar o grande tesouro final!
write_csv(tabela_final_absoluta, file.path(data_path, "Tabela_Completa_com_PROX.csv"))
message("\n >> Sucesso absoluto! A tabela com todas as métricas + PROX foi salva.")

# ==============================================================================
# SCRIPT 3: UNINDO DADOS BIOLÓGICOS, SOCIAIS, AMBIENTAIS E PAISAGEM
# ==============================================================================

library(tidyverse)
library(here)

message("--- Iniciando a Grande Fusão de Dados ---")

# 1. Carregar a tabela de métricas de paisagem (que já tem o PROX!)
tabela_paisagem <- read_csv(file.path(data_path, "Tabela_Completa_com_PROX.csv"), show_col_types = FALSE)

# 2. Carregar o seu novo arquivo com variáveis Sociais, Ambientais e Biológicas
tabela_bio_env <- read_csv(file.path(data_path, "Caracteristicas_Social_EnviromentVariables.CSV"), show_col_types = FALSE)

# 3. Preparar a tabela biológica/ambiental para o Join
tabela_bio_env_limpa <- tabela_bio_env %>%
  mutate(
    UA_ID = as.character(UA_ID),
    SPECIES = as.character(SPECIES),
    # Criamos uma chave de busca toda em minúscula para evitar erros de digitação
    SPECIES_MATCH = tolower(SPECIES) 
  )

# 4. Preparar a tabela de paisagem
tabela_paisagem_limpa <- tabela_paisagem %>%
  mutate(
    UA_ID = as.character(UA_ID),
    # Garante que temos a coluna SPECIES
    SPECIES = if("Species_Code" %in% names(.)) Species_Code else SPECIES,
    # Removemos números soltos no final do nome (ex: "mazama1" vira "mazama") e passamos pra minúscula
    SPECIES_MATCH = tolower(str_remove(SPECIES, "[0-9]+$"))
  )

# 5. O Grande Join (Unindo tudo pelo UA_ID e pela Espécie)
tabela_final_modelagem <- tabela_bio_env_limpa %>%
  inner_join(tabela_paisagem_limpa, by = c("UA_ID", "SPECIES_MATCH")) %>%
  # Limpar colunas duplicadas ou temporárias
  select(-SPECIES_MATCH, -starts_with("SPECIES.y")) %>%
  rename(SPECIES = SPECIES.x) %>% 
  # Organizar a ordem das colunas para ficar visualmente lógico
  relocate(ORDER, GENUS, SPECIES, DIET, LOCOMOTION, UA_ID, UA_Area_ha, 
           Mean_Elevation_m, Mean_Pop_Density)

# 6. Salvar a Base de Dados Definitiva!
write_csv(tabela_final_modelagem, file.path(data_path, "Data_Raw_GLM_NMDS_Final.csv"))

message("\n >> SUCESSO TOTAL! Sua base de dados 'Data_Raw_GLM_NMDS_Final.csv' está pronta para as análises estatísticas!")