# ==============================================================================
# AJUSTE FINO: GRÁFICOS PARA PUBLICAÇÃO (VERSÃO FINAL TIFF)
# ==============================================================================
library(myTAI)
library(readr)
library(ggplot2)

arquivos_limpos <- c(
  "nipponbare_full.csv", "n22_full.csv",
  "nipponbare_dusk.csv", "n22_dusk.csv",
  "nipponbare_dawn.csv", "n22_dawn.csv"
)

gerar_plot_final_v2 <- function(arquivo) {
  # 1. Identificação
  cultivar <- ifelse(grepl("nipponbare", arquivo, ignore.case = TRUE), "Nipponbare", "N22")
  cenario_nome <- gsub(".csv", "", arquivo)
  
  # 2. Carregar dados e teste
  df <- read_csv(arquivo, show_col_types = FALSE)
  phylo_set <- BulkPhyloExpressionSet_from_df(df)
  res_teste <- stat_flatline_test(phylo_set, plot_result = FALSE)
  p_val_fmt <- format.pval(res_teste@p_value, digits = 4)
  
  # 3. Construção do gráfico
  # show_p_value = FALSE remove o p-value automático do pacote
  g <- plot_signature(phylo_set, show_p_value = TRUE) +
    geom_point(size = 3.2, color = "black") + 
    labs(
      title = NULL, 
      subtitle = paste0(cultivar),
      x = "Stages",
      y = "Transcriptome age index (TAI)"
    ) +
    # Escala de Y variando em intervalos de 0.1
    scale_y_continuous(breaks = seq(0, 10, 0.1)) +
    # Estética limpa com moldura (box)
    theme_bw() + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      # Rótulos do eixo X com 45 graus
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic")
    )
  
  # 4. Salvamento em TIFF (300 DPI)
  nome_saida <- paste0("Final_Pub_", cenario_nome, ".tiff")
  ggsave(nome_saida, plot = g, width = 7, height = 5, dpi = 300, 
         device = "tiff", compression = "lzw")
  
  return(g)
}

# Execução
for(f in arquivos_limpos) {
  try(print(gerar_plot_final_v2(f)))
}
