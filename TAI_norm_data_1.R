# Librarys
library(myTAI)
library(readr)
library(ggplot2)

# Lendo os arquivos normalizados
arquivos_limpos <- c(
  "nipponbare_full_norm.csv", "n22_full_norm.csv",
  "nipponbare_dusk_norm.csv", "n22_dusk_norm.csv",
  "nipponbare_dawn_norm.csv", "n22_dawn_norm.csv"
)

gerar_plot_automatico_ajustado <- function(arquivo) {
  
  # Identificação da Cultivar
  cultivar <- ifelse(grepl("nipponbare", arquivo, ignore.case = TRUE), "Nipponbare", "N22")
  cenario_nome <- gsub(".csv", "", arquivo)
  
  message(">>> Processando: ", cenario_nome)
  
  
  #Análise 
  
  df <- read_csv(arquivo, show_col_types = FALSE)
  phylo_set <- BulkPhyloExpressionSet_from_df(df)
  res_teste <- stat_flatline_test(phylo_set, plot_result = FALSE)
  valores_tai <- TAI(phylo_set) 
  p_val_bruto <- as.numeric(res_teste@p_value)
  
  tabela_tai <- data.frame(
    Stage = names(valores_tai),
    TAI_Value = as.numeric(valores_tai),
    P_Value = p_val_bruto # Número bruto salvo na tabela
  )
  
  nome_tabela_saida <- paste0("Valores_TAI_", cenario_nome, ".csv")
  write_csv(tabela_tai, nome_tabela_saida)
  message("    -> Valores de TAI e P-valor salvos em: ", nome_tabela_saida)
  sci_val <- formatC(p_val_bruto, format = "e", digits = 2)
  sci_parts <- strsplit(sci_val, "e")[[1]]
  base <- sci_parts[1]
  exp <- as.numeric(sci_parts[2]) # Tira os zeros indesejados (ex: -05 vira -5)
  texto_p_valor <- paste0("Pflt == ", base, " %*% 10^", exp)
  
# Graph
  g <- plot_signature(phylo_set, show_p_val = FALSE)
  for (i in seq_along(g$layers)) {
    if (inherits(g$layers[[i]]$geom, "GeomRibbon") || inherits(g$layers[[i]]$geom, "GeomPolygon")) {
      g$layers[[i]]$aes_params$fill <- "#2E8B57" # Verde
    }
    if (inherits(g$layers[[i]]$geom, "GeomLine")) {
      g$layers[[i]]$aes_params$colour <- "black" # Verde
    }
  }
  g <- g + 
    geom_point(size = 3, color = "black") + 
    annotate("text", x = Inf, y = Inf, 
             label = texto_p_valor, parse = TRUE, 
             hjust = 1.2, vjust = 2, 
             size = 5, color = "black") +
    labs(
      title = NULL, 
      subtitle = cultivar, 
      x = "Stages",
      y = "Transcriptome age index (TAI)"
    ) +
    scale_y_continuous(breaks = seq(0, 10, 0.05)) +
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
  
  # Salvar a Imagem TIFF
  nome_saida <- paste0("Fig_TAI_", cenario_nome, ".tiff")
  ggsave(nome_saida, plot = g, width = 12, height = 10, dpi = 600, 
         device = "tiff", compression = "lzw")
  message("    -> Gráfico (Verde/Expoente Exato) salvo em: ", nome_saida)
  
  return(g)
}

# Executar para os 6 dataframes
for(f in arquivos_limpos) {
  try(print(gerar_plot_automatico_ajustado(f)))
}
