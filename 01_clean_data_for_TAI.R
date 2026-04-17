# ==============================================================================
# PROJETO: Filotranscritômica do Arroz (Nipponbare e N22) 
# SCRIPT 01: Limpeza de dados, médias biológicas e Dicionário de Clados.
# ==============================================================================

library(tidyverse)

# 1. Carregamento
cat("Lendo arquivos do projeto...\n")
expressao_bruta <- read_csv("norm_data_1.csv", show_col_types = FALSE)
genera_bruto    <- read_tsv("39947_gene_ages.tsv", show_col_types = FALSE)

# ==============================================================================
# 2. EXTRAÇÃO DO DICIONÁRIO DE CLADOS
# ==============================================================================
cat("Criando o Dicionário Evolutivo...\n")
dicionario_clados <- genera_bruto %>%
  # Pega só o número (rank) e o nome do clado (phylostratum)
  select(PS = rank, Clado = phylostratum) %>%
  # Mantém apenas as linhas únicas (tira as repetições dos milhares de genes)
  distinct() %>%
  # Remove NAs e organiza do 1 (mais antigo) ao 16 (mais novo)
  drop_na() %>%
  arrange(PS)

print(dicionario_clados, n = 20)

# ==============================================================================
# 3. LIMPEZA DOS IDs EVOLUTIVOS (GenEra) 
# ==============================================================================
cat("\nLimpando IDs do GenEra...\n")
genera_limpo <- genera_bruto %>%
  mutate(GeneID = ifelse(str_detect(`#gene`, "^LOC_Os"), 
                         str_replace(`#gene`, "\\..*", ""), 
                         `#gene`)) %>%
  filter(!is.na(rank)) %>%
  select(Phylostratum = rank, GeneID)

# ==============================================================================
# 4. CÁLCULO DAS MÉDIAS DAS REPLICATAS
# ==============================================================================
cat("Calculando médias...\n")
expressao_medias <- expressao_bruta %>%
  rename(GeneID = `#gene`) %>%
  pivot_longer(cols = -GeneID, names_to = c("Time", "Day", "Temp", "Genotype", "Diel", "Rep"), names_sep = "\\.", values_to = "TPM") %>%
  group_by(GeneID, Genotype, Time, Temp, Diel) %>%
  summarise(Media_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
  mutate(Estagio = paste0("T", Time, "_", Temp, "C_", Diel)) %>%
  select(GeneID, Genotype, Estagio, Media_TPM) %>%
  pivot_wider(names_from = Estagio, values_from = Media_TPM)

# ==============================================================================
# 5. INTEGRAÇÃO
# ==============================================================================
cat("Montando as matrizes finais...\n")

PhyloSet_Nipponbare <- genera_limpo %>%
  inner_join(expressao_medias %>% filter(Genotype == "Nipponbare"), by = "GeneID") %>%
  select(-Genotype)

PhyloSet_N22 <- genera_limpo %>%
  inner_join(expressao_medias %>% filter(Genotype == "N22"), by = "GeneID") %>%
  select(-Genotype)

cat("Matrizes e Dicionário gerados.\n")

# ==============================================================================
# 6. EXPORTAÇÃO DOS ARQUIVOS 
# ==============================================================================
cat("Exportando ficheiros para a pasta do projeto...\n")

# Exportar as matrizes de expressão filogenética (PhyloExpressionSets)
write_csv(PhyloSet_Nipponbare, "PhyloSet_Nipponbare.csv")
write_csv(PhyloSet_N22, "PhyloSet_N22.csv")

write_csv(expressao_medias, "Expression_cultivares_PS.csv")

# Exportar o dicionário de clados 
write_csv(dicionario_clados, "dicionario_clados.csv")

cat("Ficheiros .csv foram criados em:", getwd(), "\n")
