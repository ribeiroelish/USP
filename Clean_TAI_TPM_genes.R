# Librarys 
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# Datasets
expr_raw <- read_csv("TPM_genes.csv", show_col_types = FALSE)
meta_raw <- read_csv("Metatable.csv", show_col_types = FALSE)
gene_age_raw <- read_tsv("39947_gene_ages.tsv", show_col_types = FALSE)

# Limpeza

gene_age_clean <- gene_age_raw %>%
  rename(gene = 1, PS = 3) %>%
  filter(!is.na(PS)) %>%
  mutate(gene = if_else(str_detect(gene, "^LOC_Os"), 
                        str_remove(gene, "\\.[0-9]+$"), 
                        gene)) %>%
  group_by(gene) %>%
  summarise(PS = min(PS), .groups = "drop")

# Match PS e TPM 
meta_clean <- meta_raw %>%
  mutate(sample_id = paste(Time, Day, `Temperature.C.`, Genotype, Diel, bio_rep, sep = "."))

expr_long <- expr_raw %>%
  rename(gene = 1) %>%
  filter(!is.na(gene)) %>%
  pivot_longer(cols = -gene, names_to = "sample_id", values_to = "TPM") %>%
  inner_join(meta_clean, by = "sample_id") %>%
  inner_join(gene_age_clean, by = "gene")


# Fazer média das replicatas por stage 

data_aggregated <- expr_long %>%
  group_by(Genotype, Diel, Time, Day, gene, PS) %>%
  summarise(mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
  arrange(Time) %>%
  mutate(stage = paste0("T", Time, "_", Diel))


# Criar matriz para os 6 testes (LEMBRAR: Col1 = PS, Col2 = GeneID)


format_to_mytai <- function(df) {
  df %>%
    select(gene, PS, stage, mean_TPM) %>%
    pivot_wider(names_from = stage, values_from = mean_TPM) %>%
    rename(Phylostratum = PS, GeneID = gene) %>%
    select(Phylostratum, GeneID, everything()) %>%
    filter(rowSums(select(., -Phylostratum, -GeneID)) > 0)
}


# Todos os rstágios por cultivar
nippo_full <- data_aggregated %>% filter(Genotype == "Nipponbare") %>% format_to_mytai()
n22_full   <- data_aggregated %>% filter(Genotype == "N22") %>% format_to_mytai()

# Nipponbare (Dusk/Dawn)
nippo_dusk <- data_aggregated %>% filter(Genotype == "Nipponbare", Diel == "Dusk") %>% format_to_mytai()
nippo_dawn <- data_aggregated %>% filter(Genotype == "Nipponbare", Diel == "Dawn") %>% format_to_mytai()

# N22 (Dusk/Dawn)
n22_dusk   <- data_aggregated %>% filter(Genotype == "N22", Diel == "Dusk") %>% format_to_mytai()
n22_dawn   <- data_aggregated %>% filter(Genotype == "N22", Diel == "Dawn") %>% format_to_mytai()


# Salar 
write_csv(nippo_full, "nipponbare_full.csv")
write_csv(n22_full,   "n22_full.csv")
write_csv(nippo_dusk, "nipponbare_dusk.csv")
write_csv(nippo_dawn, "nipponbare_dawn.csv")
write_csv(n22_dusk,   "n22_dusk.csv")
write_csv(n22_dawn,   "n22_dawn.csv")


# Check de dados 
meus_arquivos <- c("nipponbare_full.csv", "n22_full.csv", 
                   "nipponbare_dusk.csv", "nipponbare_dawn.csv", 
                   "n22_dusk.csv", "n22_dawn.csv")

message("\n=== LAUDO FINAL DE NAs ===")
for(arquivo in meus_arquivos) {
  if(file.exists(arquivo)) {
    df_temp <- read.csv(arquivo)
    
    total_na <- sum(is.na(df_temp))
    
    if(total_na == 0) {
      message("[LIMPO] ", arquivo, " -> 0 NAs detectados.")
    } else {
      message("[ALERTA] ", arquivo, " -> ", total_na, " NAs encontrados!")
    }
  } else {
    message("Arquivo não encontrado: ", arquivo)
  }
}