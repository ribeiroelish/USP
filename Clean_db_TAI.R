# ==============================================================================
# 1. SETUP E CARREGAMENTO
# ==============================================================================
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# Datasets
expr_raw <- read_csv("TPM_genes.csv", show_col_types = FALSE)
meta_raw <- read_csv("Metatable.csv", show_col_types = FALSE)
gene_age_raw <- read_tsv("39947_gene_ages.tsv", show_col_types = FALSE)

# ==============================================================================
# 2. PADRONIZAÇÃO E LIMPEZA DO PHYLOSTRATUM (PS)
# ==============================================================================

gene_age_clean <- gene_age_raw %>%
  rename(gene = 1, PS = 3) %>%
  # 1º BLINDAGEM: Remove qualquer gene que não tenha a idade anotada (NA)
  filter(!is.na(PS)) %>%
  # 2º BLINDAGEM: Limpeza dos sufixos (.1, .2) APENAS para os genes LOC_Os
  mutate(gene = if_else(str_detect(gene, "^LOC_Os"), 
                        str_remove(gene, "\\.[0-9]+$"), 
                        gene)) %>%
  # Agrupa as isoformas limpas e mantém a idade mais antiga (menor PS)
  group_by(gene) %>%
  summarise(PS = min(PS), .groups = "drop")

# ==============================================================================
# 3. PREPARAÇÃO DA EXPRESSÃO E O CRUZAMENTO (MATCH)
# ==============================================================================
meta_clean <- meta_raw %>%
  mutate(sample_id = paste(Time, Day, `Temperature.C.`, Genotype, Diel, bio_rep, sep = "."))

expr_long <- expr_raw %>%
  rename(gene = 1) %>%
  # 3º BLINDAGEM: Remove possíveis NAs no nome do gene na tabela de expressão
  filter(!is.na(gene)) %>%
  pivot_longer(cols = -gene, names_to = "sample_id", values_to = "TPM") %>%
  inner_join(meta_clean, by = "sample_id") %>%
  # O inner_join ABAIXO garante que SÓ os genes correspondentes em ambas as tabelas sobrevivam
  inner_join(gene_age_clean, by = "gene")

# ==============================================================================
# 4. MÉDIA DAS RÉPLICAS BIOLÓGICAS
# ==============================================================================
data_aggregated <- expr_long %>%
  group_by(Genotype, Diel, Time, Day, gene, PS) %>%
  summarise(mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
  arrange(Time) %>%
  mutate(stage = paste0("T", Time, "_", Diel))

# ==============================================================================
# 5. MATRIZ 6 CENÁRIOS PARA myTAI (FORMATO S7 ESTRITO)
# ==============================================================================

format_to_mytai <- function(df) {
  df %>%
    select(gene, PS, stage, mean_TPM) %>%
    pivot_wider(names_from = stage, values_from = mean_TPM) %>%
    rename(Phylostratum = PS, GeneID = gene) %>%
    # Garante a ordem exata exigida pelo pacote: Col1 = PS, Col2 = GeneID
    select(Phylostratum, GeneID, everything()) %>%
    # 4º BLINDAGEM: Remove genes que, após a média, ficaram com expressão total zero
    filter(rowSums(select(., -Phylostratum, -GeneID)) > 0)
}

# Criando os subsets
# 1 & 2: Full (Todos os pontos por cultivar)
nippo_full <- data_aggregated %>% filter(Genotype == "Nipponbare") %>% format_to_mytai()
n22_full   <- data_aggregated %>% filter(Genotype == "N22") %>% format_to_mytai()

# 3 & 4: Nipponbare separado por Diel
nippo_dusk <- data_aggregated %>% filter(Genotype == "Nipponbare", Diel == "Dusk") %>% format_to_mytai()
nippo_dawn <- data_aggregated %>% filter(Genotype == "Nipponbare", Diel == "Dawn") %>% format_to_mytai()

# 5 & 6: N22 separado por Diel
n22_dusk   <- data_aggregated %>% filter(Genotype == "N22", Diel == "Dusk") %>% format_to_mytai()
n22_dawn   <- data_aggregated %>% filter(Genotype == "N22", Diel == "Dawn") %>% format_to_mytai()


# ==============================================================================
# 6. EXPORTAÇÃO
# ==============================================================================
write_csv(nippo_full, "nipponbare_full.csv")
write_csv(n22_full,   "n22_full.csv")
write_csv(nippo_dusk, "nipponbare_dusk.csv")
write_csv(nippo_dawn, "nipponbare_dawn.csv")
write_csv(n22_dusk,   "n22_dusk.csv")
write_csv(n22_dawn,   "n22_dawn.csv")


# Lista dos arquivos que você gerou
meus_arquivos <- c("nipponbare_full.csv", "n22_full.csv", 
                   "nipponbare_dusk.csv", "nipponbare_dawn.csv", 
                   "n22_dusk.csv", "n22_dawn.csv")

message("\n=== LAUDO FINAL DE NAs ===")
for(arquivo in meus_arquivos) {
  if(file.exists(arquivo)) {
    # Lemos como dataframe puro para varrer tudo
    df_temp <- read.csv(arquivo)
    
    # Soma de NAs em todas as linhas e colunas
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