# ==============================================================================
# PROJETO: Filotranscritômica do Arroz - VERSÃO COM REPLICATAS
# ==============================================================================

library(tidyverse)

# 1. Carregamento
cat("Lendo arquivos brutos...\n")
expressao_bruta <- read_csv("norm_data_1.csv", show_col_types = FALSE)
genera_bruto    <- read_tsv("39947_gene_ages.tsv", show_col_types = FALSE)

# 2. Dicionário de Clados 
dicionario_clados <- genera_bruto %>%
  select(PS = rank, Clado = phylostratum) %>%
  distinct() %>%
  drop_na() %>%
  arrange(PS)

# 3. Limpeza dos IDs (GenEra)
genera_limpo <- genera_bruto %>%
  mutate(GeneID = ifelse(str_detect(`#gene`, "^LOC_Os"), 
                         str_replace(`#gene`, "\\..*", ""), 
                         `#gene`)) %>%
  filter(!is.na(rank)) %>%
  select(Phylostratum = rank, GeneID)

# ==============================================================================
# 4. ORGANIZAÇÃO DAS REPLICATAS 
# ==============================================================================
cat("Organizando as replicatas biológicas...\n")

expressao_longa <- expressao_bruta %>%
  rename(GeneID = `#gene`) %>%
  pivot_longer(
    cols = -GeneID,
    names_to = c("Time", "Day", "Temp", "Genotype", "Diel", "Rep"),
    names_sep = "\\.",
    values_to = "TPM"
  ) %>%

  mutate(Estagio_Rep = paste0("T", Time, "_", Temp, "C_", Diel, "_", Rep))

# ==============================================================================
# 5. MONTAGEM DAS MATRIZES (PhyloExpressionSets)
# ==============================================================================
cat("Gerando matrizes para Nipponbare e N22...\n")

# Matriz Nipponbare
PhyloSet_Nipponbare_Reps <- genera_limpo %>%
  inner_join(
    expressao_longa %>% 
      filter(Genotype == "Nipponbare") %>%
      select(GeneID, Estagio_Rep, TPM) %>%
      pivot_wider(names_from = Estagio_Rep, values_from = TPM),
    by = "GeneID"
  )

# Matriz N22
PhyloSet_N22_Reps <- genera_limpo %>%
  inner_join(
    expressao_longa %>% 
      filter(Genotype == "N22") %>%
      select(GeneID, Estagio_Rep, TPM) %>%
      pivot_wider(names_from = Estagio_Rep, values_from = TPM),
    by = "GeneID"
  )

# ==============================================================================
# 6. EXPORTAÇÃO
# ==============================================================================
write_csv(PhyloSet_Nipponbare_Reps, "PhyloSet_Nipponbare_Reps.csv")
write_csv(PhyloSet_N22_Reps, "PhyloSet_N22_Reps.csv")
write_csv(dicionario_clados, "dicionario_clados.csv")


