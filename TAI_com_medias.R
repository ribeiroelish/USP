# ==============================================================================
# SCRIPT 02: Análise Filotranscritômica TAI (Nipponbare e N22)
# ==============================================================================

library(tidyverse)
library(myTAI)

# 1. Criar função 
processar_tai <- function(arquivo_csv, nome_variedade) {
  
  cat("Processando:", nome_variedade, "...\n")
  
  # A. Carregar e transformar em objeto myTAI
  df <- read_csv(arquivo_csv, show_col_types = FALSE)
  phylo_set <- BulkPhyloExpressionSet_from_df(df)
  
  # B. Gerar o gráfico
  p <- plot_signature(
    phyex_set  = phylo_set,
    show_reps  = FALSE,
    show_p_val = FALSE
  ) + 
    labs(
      title = paste(nome_variedade),
      x = "Stages",
      y = "Transcriptome Age Index (TAI)"
    ) +
    theme_minimal()
  
  # C. Salvar a imagem 
  nome_arquivo <- paste0("Plot_TAI_Medias_", nome_variedade, ".png")
  ggsave(nome_arquivo, plot = p, width = 8, height = 5, dpi = 300)
  
  cat("Gráfico salvo como:", nome_arquivo, "\n\n")
  
  return(p)
}

# Por cultivar 
plot_nipponbare <- processar_tai("PhyloSet_Nipponbare.csv", "Nipponbare")
plot_n22         <- processar_tai("PhyloSet_N22.csv", "N22")

print(plot_nipponbare)
print(plot_n22)
