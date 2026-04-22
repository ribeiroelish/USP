# Librarys

library(myTAI)
library(readr)
library(ggplot2)

arquivos_limpos <- c(
  "nipponbare_full.csv", "n22_full.csv",
  "nipponbare_dusk.csv", "n22_dusk.csv",
  "nipponbare_dawn.csv", "n22_dawn.csv"
)

gerar_plot_automatico_ajustado <- function(arquivo) {
 
# Identificação da Cultivar
  
  cultivar <- ifelse(grepl("nipponbare", arquivo, ignore.case = TRUE), "Nipponbare", "N22")
  cenario_nome <- gsub(".csv", "", arquivo)
  
  message(">>> Processando: ", cenario_nome)
  
  
  df <- read_csv(arquivo, show_col_types = FALSE)
  phylo_set <- BulkPhyloExpressionSet_from_df(df)
  res_teste <- stat_flatline_test(phylo_set, plot_result = FALSE)
  
  texto_p_valor <- paste0("Pflt = ", format.pval(res_teste@p_value, digits = 3))
  
  
  g <- plot_signature(phylo_set, show_p_val = FALSE) +
    geom_point(size = 3, color = "black") + 
    
    #p-value anotado 
    annotate("text", x = Inf, y = Inf, 
             label = texto_p_valor, 
             hjust = 1.2, vjust = 2, 
             size = 5, color = "black") +
    labs(
      title = NULL, 
      subtitle = cultivar, 
      x = "Stages",
      y = "Transcriptome age index (TAI)"
    ) +
    scale_y_continuous(breaks = seq(0, 10, 0.1)) +
    theme_bw() + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic")
    )
  
  # Salvar
  nome_saida <- paste0("Final_Automatico_", cenario_nome, ".tiff")
  ggsave(nome_saida, plot = g, width = 12, height = 10, dpi = 600, 
         device = "tiff", compression = "lzw")
  
  return(g)
}

# Executar pra os 6 df
for(f in arquivos_limpos) {
  try(print(gerar_plot_automatico_ajustado(f)))
}
