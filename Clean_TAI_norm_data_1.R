
# Librarys 
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# Datasets 
expr_raw <- read_csv("norm_data_1.csv", show_col_types = FALSE)
gene_age_raw <- read_tsv("39947_gene_ages.tsv", show_col_types = FALSE)

# Clean
gene_age_clean <- gene_age_raw %>%
  rename(gene = 1, PS = 3) %>%
  filter(!is.na(PS)) %>%
  mutate(gene = if_else(str_detect(gene, "^LOC_Os"), 
                        str_remove(gene, "\\.[0-9]+$"), 
                        gene)) %>%
  group_by(gene) %>%
  summarise(PS = min(PS), .groups = "drop")

# Match PS e Expressão 
expr_long <- expr_raw %>%
  rename(gene = 1) %>%
  mutate(gene = if_else(str_detect(gene, "^LOC_Os"), 
                        str_remove(gene, "\\.[0-9]+$"), 
                        gene)) %>%
  filter(!is.na(gene)) %>%
  pivot_longer(cols = -gene, names_to = "sample_id", values_to = "norm_expr") %>%
  separate(sample_id, 
           into = c("Time", "Day", "Temp", "Genotype", "Diel", "Rep"), 
           sep = "\\.") %>%
  mutate(Time = as.numeric(Time)) %>%
  inner_join(gene_age_clean, by = "gene")

# Fazer média das replicatas por stage 
data_aggregated <- expr_long %>%
  group_by(Genotype, Diel, Time, Day, gene, PS) %>%
  summarise(mean_expr = mean(norm_expr, na.rm = TRUE), .groups = "drop") %>%
  arrange(Time) %>%
  mutate(stage = paste0("T", Time, "_", Diel))

# Função para formatar as matrizes para os 6 testes 
format_to_mytai <- function(df) {
  df %>%
    select(gene, PS, stage, mean_expr) %>%
    pivot_wider(names_from = stage, values_from = mean_expr) %>%
    rename(Phylostratum = PS, GeneID = gene) %>%
    select(Phylostratum, GeneID, everything()) %>%
    filter(rowSums(select(., -Phylostratum, -GeneID)) > 0)
}

# Todos os estágios por cultivar
nippo_full <- data_aggregated %>% filter(Genotype == "Nipponbare") %>% format_to_mytai()
n22_full   <- data_aggregated %>% filter(Genotype == "N22") %>% format_to_mytai()

# Nipponbare (Dusk/Dawn)
nippo_dusk <- data_aggregated %>% filter(Genotype == "Nipponbare", Diel == "Dusk") %>% format_to_mytai()
nippo_dawn <- data_aggregated %>% filter(Genotype == "Nipponbare", Diel == "Dawn") %>% format_to_mytai()

# N22 (Dusk/Dawn)
n22_dusk   <- data_aggregated %>% filter(Genotype == "N22", Diel == "Dusk") %>% format_to_mytai()
n22_dawn   <- data_aggregated %>% filter(Genotype == "N22", Diel == "Dawn") %>% format_to_mytai()


# Salvar 
write_csv(nippo_full, "nipponbare_full_norm.csv")
write_csv(n22_full,   "n22_full_norm.csv")
write_csv(nippo_dusk, "nipponbare_dusk_norm.csv")
write_csv(nippo_dawn, "nipponbare_dawn_norm.csv")
write_csv(n22_dusk,   "n22_dusk_norm.csv")
write_csv(n22_dawn,   "n22_dawn_norm.csv")


# Check de dados finais
meus_arquivos <- c("nipponbare_full_norm.csv", "n22_full_norm.csv", 
                   "nipponbare_dusk_norm.csv", "nipponbare_dawn_norm.csv", 
                   "n22_dusk_norm.csv", "n22_dawn_norm.csv")

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
