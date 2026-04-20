# ==============================================================================
# PROJETO: Filotranscritômica - PREPARAÇÃO DE DADOS (Dusk/Down)
# ==============================================================================

library(tidyverse)

# 1. Datasets
expressao_bruta <- read_csv("TPM_genes.csv", show_col_types = FALSE)
genera_bruto    <- read_tsv("39947_gene_ages.tsv", show_col_types = FALSE)

# 2. Dicionário de Clados 
dicionario_clados <- genera_bruto %>%
  select(PS = rank, Clado = phylostratum) %>%
  distinct() %>%
  drop_na() %>%
  arrange(PS)

# 3. Clean

genera_limpo <- genera_bruto %>%
  mutate(GeneID = ifelse(
    str_detect(`#gene`, "^LOC_Os"), 
    str_replace(`#gene`, "\\.\\d+$", ""), 
    `#gene`
  )) %>%
  filter(!is.na(rank)) %>%
  distinct(GeneID, .keep_all = TRUE) %>%
  select(Phylostratum = rank, GeneID)

# ==============================================================================
# 4. ORGANIZAÇÃO DAS REPLICATAS 
# ==============================================================================

expressao_longa <- expressao_bruta %>%
  rename(GeneID = `#gene`) %>%
  pivot_longer(
    cols = -GeneID,
    names_to = c("Time", "Day", "Temp", "Genotype", "Diel", "Rep"),
    names_sep = "\\.",
    values_to = "TPM"
  ) %>%
  mutate(Estagio_Rep = paste0("T", Time, "_", Temp, "C_", Rep))

# ==============================================================================
# 5.(Nipponbare/N22 x Dusk/Dawn)
# ==============================================================================


# Função 
criar_matriz_tai <- function(df_longo, idades, genotipo, turno) {
  idades %>%
    inner_join(
      df_longo %>% 
        filter(Genotype == genotipo, Diel == turno) %>%
        arrange(as.numeric(Time), Rep) %>%
        select(GeneID, Estagio_Rep, TPM) %>%
        pivot_wider(names_from = Estagio_Rep, values_from = TPM),
      by = "GeneID"
    )
}

# Quatro grupos
PhyloSet_Nip_Dusk   <- criar_matriz_tai(expressao_longa, genera_limpo, "Nipponbare", "Dusk")
PhyloSet_Nip_Dawn <- criar_matriz_tai(expressao_longa, genera_limpo, "Nipponbare", "Dawn")

PhyloSet_N22_Dusk   <- criar_matriz_tai(expressao_longa, genera_limpo, "N22", "Dusk")
PhyloSet_N22_Dawn <- criar_matriz_tai(expressao_longa, genera_limpo, "N22", "Dawn")

# ==============================================================================
# 6. EXPORTAÇÃO
# ==============================================================================

write_csv(PhyloSet_Nip_Dusk, "PhyloSet_Nip_Dusk.csv")
write_csv(PhyloSet_Nip_Dawn, "PhyloSet_Nip_Dawn.csv")
write_csv(PhyloSet_N22_Dusk, "PhyloSet_N22_Dusk.csv")
write_csv(PhyloSet_N22_Dawn, "PhyloSet_N22_Dawn.csv")
write_csv(dicionario_clados, "dicionario_clados.csv")




