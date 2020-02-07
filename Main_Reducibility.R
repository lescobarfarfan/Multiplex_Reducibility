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

path_to_write_output <- "Write the path where output should be stores:\n" %>%
    readline()
path_to_multiplex <- "Write the path to the .Rds file containing the multiplex:\n" %>%
    readline()

# path_current_file <- rstudioapi::getActiveDocumentContext()$path
# path_current_file <- str_split(string = path_current_file, pattern = "/")
# path_current_file <- unlist(path_current_file)
# path_current_file <- path_current_file[1:(length(path_current_file) - 1)]
# path_current_file <- paste0(path_current_file, collapse = "/")
# path_current_file <- paste0(path_current_file, "/")

path_to_aux_functions <- "Write the path to the directory with the auxiliary functions" %>%
    readline()


# Source script containg necessary functions ------------------------------
# source_files_path <- paste0(path_to_aux_functions, "AuxiliaryFunctions/")
source.files <- list.files(path_to_aux_functions)
map(.x = paste0(path_to_aux_functions, source.files), .f = source, verbose = F)

# Load multiplex network previously created
multicapa_temporal <- readRDS(path_to_multiplex)


# Number of layers
no_layers <- length(multicapa_temporal[[1]]) - 1
names_layers <- names(multicapa_temporal[[1]])
names_layers <- names_layers[-length(names_layers)]

# Perform optimization for each date --------------------------------------


fechas <- names(multicapa_temporal)
fechas <- as.Date(fechas)

tic()
optimal_temporal <- pmap(.l = list(multilayer = multicapa_temporal,
                                   name_agg = "Agregada", directed = T,
                                   weighted = T),
                         .f = greedy_reduction)
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
    scale_x_date(date_breaks = "6 month", date_labels = format("%Y-%m"))


ggsave(filename = "HistoricalReducibility.pdf", plot = red_plot, device = "pdf",
       path = path_to_write_output, scale = 2, width = 19, height = 9.5,
       units = "cm", dpi = 300)



#   -----------------------------------------------------------------------
# Calculate structural metrics of multiplex -------------------------------
#   -----------------------------------------------------------------------


# Multiplex degree --------------------------------------------------------
grados_temporales <- map(.x = multicapa_temporal,
                         .f = calculate_multiplex_degree,
                         directed = F)


grados_temporales_df <- grados_temporales %>%
    bind_rows(.id = "Fecha")
grados_temporales_df$Fecha <- grados_temporales_df$Fecha %>% as.Date()

plot_degrees_temporal_by_layer <- grados_temporales_df %>%
    gather(key = Capa, value = Grado, -Fecha, -Node) %>%
    filter(Capa != "Agregada",
           Capa != "OverlappingDeg") %>%
    group_by(Fecha, Capa) %>%
    summarize(`Minimum degree` = min(Grado, na.rm = T),
              `Mean degree` = mean(Grado, na.rm = T),
              `Maximun degree` = max(Grado, na.rm = T)) %>%
    ungroup() %>%
    gather(key = Serie, value = Valor, -Fecha, -Capa) %>%
    ggplot(aes(x = Fecha, y = Valor, group = Serie, colour = Serie)) +
    geom_ma(ma_fun = SMA, n = 20, linetype = "solid") +
    geom_vline(xintercept = as.Date("2018-07-02"), linetype = "longdash",
               colour = "red", alpha = 0.7) +
    theme_bw() +
    theme(axis.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 13),
          legend.title = element_blank(),
          axis.text = element_text(size = 13),
          axis.text.x = element_text(angle = 90),
          strip.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 14)) +
    facet_wrap(~Capa) +
    scale_colour_viridis_d(option = "D") +
    scale_x_date(date_breaks = "6 months", date_labels = format("%Y-%m")) +
    labs(title = "20-day simple moving average of max, min and mean degree per layer",
         subtitle = "Number of connections")

ggsave(filename = "HistoricalDegreeByLayer.pdf",
       plot = plot_degrees_temporal_by_layer, device = "pdf",
       path = path_to_write_output, scale = 2, width = 19, height = 9.5,
       units = "cm", dpi = 300)







# Node activity -----------------------------------------------------------
node_activity <- grados_temporales_df %>%
    select(-Agregada, -OverlappingDeg) %>%
    gather(key = Capa, value = Valor, -Fecha, -Node) %>%
    mutate(Activo = if_else(Valor > 0, 1, 0))

node_activity <- node_activity %>%
    select(-Valor) %>%
    spread(key = Capa, value = Activo)



# Active nodes per layer --------------------------------------------------
actives_nodes_per_layer <- node_activity %>%
    gather(key = Capa, value = Status, -Fecha, -Node) %>%
    group_by(Fecha, Capa) %>%
    summarize(NoActiveNodes = sum(Status, na.rm = T)) %>%
    ggplot(aes(x = Fecha, y = NoActiveNodes, group = Capa, colour = Capa)) +
    geom_ma(ma_fun = SMA, n = 20, linetype = "solid") +
    theme_bw() +
    theme(axis.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 13),
          legend.title = element_blank(),
          axis.text = element_text(size = 13),
          axis.text.x = element_text(angle = 90),
          strip.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 14)) +
    scale_color_viridis_d(option = "D") +
    scale_x_date(date_breaks = "6 months", date_labels = format("%Y-%m")) +
    labs(title = "20-days moving average of the number of active banks per layer",
         subtitle = "Number of institutions")

ggsave(filename = "NoActivesNodesPerLayer.pdf", plot = actives_nodes_per_layer,
       device = "pdf", path = path_to_write_output, scale = 2, width = 19,
       height = 9.5, units = "cm", dpi = 400)



# Active layers per node --------------------------------------------------

active_layers_per_node <- node_activity %>%
    gather(key = Capa, value = Status, -Fecha, -Node) %>%
    group_by(Fecha, Node) %>%
    summarize(NoActiveLayers = sum(Status, na.rm = T)) %>%
    ggplot(aes(x = Fecha, y = NoActiveLayers, group = Node, colour = Node)) +
    geom_ma(ma_fun = SMA, n = 20, linetype = "solid") +
    geom_vline(xintercept = as.Date("2018-07-02"), colour = "red",
               linetype = "longdash") +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.title = element_blank(),
          axis.text = element_text(size = 13),
          axis.text.x = element_text(angle = 90),
          strip.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 14)) +
    facet_wrap(~Node, scales = "free_y") +
    scale_x_date(date_breaks = "6 months", date_labels = format("%Y-%m")) +
    scale_color_viridis_d(option = "D") +
    labs(title = "20-days moving average of the number of active layers per node",
         subtitle = "Number of markets")

ggsave(filename = "NoActiveLayersPerNode.pdf", plot = active_layers_per_node,
       device = "pdf", path = path_to_write_output, scale = 2, width = 19,
       height = 9.5, units = "cm", dpi = 400)





# Layer participation for each node ---------------------------------------
layer_participation_coefficient_per_node <- grados_temporales_df %>%
    select(-Agregada) %>%
    mutate_at(vars(-c(Fecha, Node, OverlappingDeg)),
              funs(./OverlappingDeg))


layer_participation_coefficient_per_node <- layer_participation_coefficient_per_node %>%
    select(-OverlappingDeg) %>%
    mutate_at(vars(-c(Fecha, Node)), funs(.^(2)))


layer_participation_coefficient_per_node <- layer_participation_coefficient_per_node %>%
    mutate(SumaAux = rowSums(select(layer_participation_coefficient_per_node, -Fecha, -Node)))
