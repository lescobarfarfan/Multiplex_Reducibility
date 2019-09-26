#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
# Main script to calculate reducibility matrices --------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(Matrix)
library(tictoc)
library(tidyquant)



# Paths -------------------------------------------------------------------

path_to_write_output <- "C:/Users/K15523/Documents/Trabajo/Multiplex_Reducibility/Resultados/"

# Source script containg necessary functions ------------------------------
source_files_path <- "C:/Users/K15523/Desktop/Multiplex_Reducibility/AuxiliaryFunctions/"
source.files <- list.files(source_files_path)
map(.x = paste0(source_files_path, source.files), .f = source)

# Load multiplex network previously created
multicapa_temporal <- readRDS("C:/Users/K15523/Documents/Trabajo/Multiplex_Reducibility/Resultados/MulticapaTemporal_Ad.Rds")



# Perform optimization for each date --------------------------------------


fechas <- names(multicapa_temporal)

tic()
optimal_temporal <- map(.x = multicapa_temporal, .f = greedy_reduction)
toc()

# Save the output
saveRDS(object = optimal_temporal,
        file = paste0(path_to_write_output, "Optimal_Aggregations.Rds"))


# Obtain a plot of the historical reducibility ----------------------------
reducibilidad_hist <- list()
for (i in 1:length(optimal_temporal)) {
    reducibilidad_hist[i] <- optimal_temporal[[i]]$Reducibility
}
names(reducibilidad_hist) <- fechas

reducibilidad_df <- reducibilidad_hist %>% bind_rows()
reducibilidad_df <- reducibilidad_df %>%
    gather(key = Fecha, value = Red)


reducibilidad_df$Fecha <- reducibilidad_df$Fecha %>% as.Date()
fechas <- fechas %>% as.Date()

red_plot <- reducibilidad_df %>%
    ggplot(aes(x = Fecha, y = Red, colour = "red")) +
    geom_ma(ma_fun = SMA, n = 21, linetype = "solid") +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 18, face = "bold"),
          plot.subtitle = element_text(size = 15),
          axis.text.x = element_text(angle = 90)) +
    ggtitle(label = "20-day simple moving average of the reducibility of the system",
            subtitle = paste0("From ", fechas[1], " to ", fechas[2])) +
    scale_x_date(date_breaks = "1 month", date_labels = format("%Y-%m"))


ggsave(filename = "HistoricalReducibility.pdf", plot = red_plot, device = "pdf",
       path = path_to_write_output, scale = 2, width = 19, height = 9.5,
       units = "cm", dpi = 300)
