# ==============================================================================
# SCRIPT 02: Análise Filotranscritômica TAI (Versão Estatística)
# OBJETIVO: Calcular TAI com variância e Teste de Significância (FlatLineTest)
# ==============================================================================

library(tidyverse)
library(myTAI)

# 1. FUNÇÃO PARA PROCESSAMENTO ESTATÍSTICO
processar_tai_estatistico <- function(arquivo_csv, nome_variedade) {
  
  cat("Lendo dados e calculando variância para:", nome_variedade, "...\n")
  
  # A. Carregar a matriz com replicatas
  df_reps <- read_csv(arquivo_csv, show_col_types = FALSE)
  
  # B. Converter para o objeto moderno (S7) do myTAI
  # O pacote agrupa automaticamente as colunas que pertencem ao mesmo estágio
  phylo_reps <- BulkPhyloExpressionSet_from_df(df_reps)
  
  # C. Gerar o gráfico com Replicatas e P-value
  # show_reps = TRUE desenha a dispersão das triplicatas
  # show_p_val = TRUE executa o FlatLineTest por padrão
  p <- plot_signature(
    phyex_set  = phylo_reps,
    show_reps  = TRUE,      
    show_p_val = TRUE       
  ) + 
    labs(
      title = paste(nome_variedade),
      x = "Stages",
      y = "Transcriptome Age Index (TAI)"
    ) +
    theme_minimal()
  

  nome_img <- paste0("Plot_TAI_ESTATISTICO_", nome_variedade, ".png")
  ggsave(nome_img, plot = p, width = 10, height = 6, dpi = 300)
  
  
  return(p)
}

# 2. EXECUÇÃO PARA NIPPONBARE E N22
# Usando os arquivos gerados com as triplicatas
plot_nippon_reps <- processar_tai_estatistico("PhyloSet_Nipponbare_Reps.csv", "Nipponbare")
plot_n22_reps     <- processar_tai_estatistico("PhyloSet_N22_Reps.csv", "N22")

print(plot_nippon_reps)
print(plot_n22_reps)

