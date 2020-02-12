#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
# Script to calculate multilayer debtrank as in Poledna et al -------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------


require(NetworkRiskMeasures)
require(igraph)
require(tidyverse)
require(tidyquant)

# Load the data
path_to_multilayer <- "/media/luisescobar/TRABAJO/Multilayer/Results/MulticapaTemporal_Contagio.Rds"
path_to_mult_reduced <- "/media/luisescobar/TRABAJO/Multilayer/Results/Optimal_Aggregations.Rds"
path_to_thresholds <- "/media/luisescobar/TRABAJO/Multilayer/Umbrales/"

# Load the multiplex network
multiplex <- readRDS(path_to_multilayer)

dates <- names(multiplex)
dates <- dates %>% str_sub(start = 3, end = 10)
dates <- dates %>% str_remove_all(pattern = "-")


# Load the assets and the capital of banks
files <- list.files(path = path_to_thresholds, pattern = "*.csv")

# Intersection of files and dates in multiplex
files_2 <- files %>% str_sub(start = 1, end = 6)
intersection <- dplyr::intersect(dates, files_2)

files_3 <- files[which(str_sub(files, 1, 6) %in% intersection)]

# Load the assets and capital data
extra_data <- map(.x = paste0(path_to_thresholds, files_3), 
                  .f = read_csv,
                  col_names = F)
names(extra_data) <- names(multiplex)


for (i in 1:length(extra_data)) {
  names(extra_data[[i]]) <- c("Capital", "RWA", "Assets", "Threshold")
}


#' Function to calculate systemic risk using the DebtRank following the
#' methodology presented in Poledna et. al. (2015). This function takes advantage
#' of the contagion function from NetworkRiskMeasures package.
#' @param multicapa A list containing all the layers from a multilayer network.
#' @param assets A vector containing the total assets for each node in the
#' multiplex
#' @param capital A vector containing the capital buffer for each node in the
#' multiplex
#' @param normalized If TRUE, the function will return the normalized DebtRank
#' por each layer proposed in Poledna et. al. (2015).
#' @param ... aditional parameters for the contagion function (see contagion
#' function from NetworkRiskMeasures package)
#' @return A list containing the contribution to systemic risk for each layer,
#' as in Poledna et. al. (2015)
#' @import NetworkRiskMeasures
#' @import igraph
#' @import purrr
#' @export
multiplex_systemic_risk <- function(multilayer, assets, capital,
                                    normalized = TRUE) {
  if (any(assets == 0)) {
    assets[assets == 0] <- 0.0001
  }
  if (any(capital == 0)) {
    capital[capital == 0] <- 0.0001
  }
  # First, we calculate the debt rank for each node for each layer
  
  V <- map(.x = multilayer, .f = sum)
  R <- list()
  for (i in 1:length(multilayer)) {
    layer <- names(multilayer)[i]
    aux_dr <- NetworkRiskMeasures::contagion(exposures = multilayer[[i]],
                                             buffer = capital, weights = assets,
                                             shock = "all", method = "debtrank",
                                             verbose = F)
    aux_dr_2 <- summary(aux_dr)$summary_table %>% tbl_df()
    aux_dr_2 <- aux_dr_2 %>% mutate(DebtRank = original_stress + additional_stress)
    aux_dr_2 <- aux_dr_2 %>% select(scenario, DebtRank) %>% rename(Institution = scenario)
    aux_dr_2$Institution <- aux_dr_2$Institution %>% as.character() %>% as.integer()
    aux <- aux_dr_2$DebtRank
    names(aux) <- aux_dr_2$Institution
    R[[layer]] <- aux
  }
  # Now, in order to make comparable the debtrank across layers, we use the
  # normalization proposed in Poledna et. al. (2015).
  R_norm <- list()
  for (i in 1:length(R)) {
    layer <- names(R)[i]
    aux_v <- V[[layer]]/V[["Agregada"]]
    R_norm[[layer]] <- R[[layer]]*aux_v
  }
  if (normalized) {
    return(R_norm)
  } else {
    return(R)
  }
}



# Analysis
resultados <- data.frame() %>% tbl_df()
for (i in 1:length(multiplex)) {
  cat("Processing date ", i, " of ", length(multiplex), "\n\n")
  aux_res <- multiplex_systemic_risk(multilayer = multiplex[[i]],
                                     assets = extra_data[[i]]$Assets,
                                     capital = extra_data[[i]]$Capital,
                                     normalized = T)
  institutions <- names(aux_res[[1]])
  aux_res2 <- aux_res %>% bind_cols()
  aux_res2$Institution <- institutions
  aux_res2 <- aux_res2 %>% 
    select(Institution, 1:(ncol(aux_res2) - 1))
  aux_res2$Date = names(multiplex)[i]
  aux_res2 <- aux_res2 %>%
    select(Date, 1:(ncol(aux_res2) - 1))
  aux_res2$Date <- aux_res2$Date %>% as.Date()
  resultados <- list(resultados, aux_res2) %>% bind_rows()
}
# Multiply by a 100 to make values percentages
resultados_gather <- resultados %>% 
  gather(key = Layer, value = DebtRank, -Date, -Institution)
resultados_gather <- resultados_gather %>%
  mutate(DebtRank = DebtRank*100)


# Save results
resultados_gather %>% write_csv(path = "/media/luisescobar/TRABAJO/Multilayer/Results/DebtRankPerLayer.csv",
                                col_names = T)



# Average DebtRank per layer ----------------------------------------------
avg_dr_per_layer <- resultados_gather %>%
  group_by(Date, Layer) %>%
  summarize(AvgDebtRank = mean(DebtRank, na.rm = T)) %>%
  ungroup()

avg_dr_per_layer_minus_agg <- avg_dr_per_layer %>%
  filter(Layer != "Agregada")

avg_dr_per_layer_agg <- avg_dr_per_layer %>%
  filter(Layer == "Agregada")


avg_dr_per_layer_minus_agg %>%
  ggplot(aes(x = Date, y = AvgDebtRank, fill = Layer)) +
  geom_area() +
  geom_line(mapping = aes(x = Date, y = AvgDebtRank, colour = Layer),
            data = avg_dr_per_layer_agg) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_breaks = "1 year", date_labels = format("%Y-%m")) +
  labs(title = "Normalized DebtRank per layer in the multiplex exposures network",
       subtitle = "21-days simple moving average")
  


# Maximum DebtRank per layer ----------------------------------------------
# resultados_gather %>%
#   group_by(Date, Layer) %>%
#   summarize(MaxDebtRank = max(DebtRank, na.rm = T)) %>%
#   ggplot(aes(x = Date, y = MaxDebtRank, colour = Layer)) +
#   # geom_ma(ma_fun = SMA, n = 21, linetype = "solid")+
#   geom_line() +
#   theme_bw() +
#   theme(axis.title = element_blank(),
#         axis.text = element_text(size = 14),
#         legend.position = "bottom",
#         legend.title = element_blank(),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, face = "bold"),
#         plot.subtitle = element_text(size = 14),
#         axis.text.x = element_text(angle = 90)) +
#   scale_x_date(date_breaks = "1 year", date_labels = format("%Y-%m")) +
#   labs(title = "Maximum normalized DebtRank per layer in the multiplex exposures network",
#        subtitle = "21-days simple moving average")





# Average DebtRank per institution ----------------------------------------
avg_dr_per_inst <- resultados_gather %>%
  group_by(Institution, Layer) %>%
  summarize(MeanDebtRank = mean(DebtRank, na.rm = T))

avg_dr_per_inst_minus_agg <- avg_dr_per_inst %>%
  filter(Layer != "Agregada")


avg_dr_per_inst_agg <- avg_dr_per_inst %>%
  filter(Layer == "Agregada") %>%
  arrange(desc(MeanDebtRank))

order_inst <- avg_dr_per_inst_agg$Institution

avg_dr_per_inst_agg$Institution <- avg_dr_per_inst_agg$Institution %>% factor(levels = order_inst,
                                                                              ordered = T)
avg_dr_per_inst_minus_agg$Institution <- avg_dr_per_inst_minus_agg$Institution %>%
  factor(levels = order_inst, ordered = T)

avg_dr_per_inst_minus_agg %>%
  ggplot(aes(x = Institution, y = MeanDebtRank, fill = Layer)) +
  geom_col() +
  geom_step(aes(x = Institution, y = MeanDebtRank, group=1, colour = Layer), data = avg_dr_per_inst_agg) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = c(0.85, 0.75),
        panel.grid = element_blank()) +
  labs(title = "Mean DebtRank per institution across the entire period of analysis",
       subtitle = "Percentage of assets of the system under stress")







# Analysis with the reduced multiplex system ------------------------------
multiplex_red <- readRDS(path_to_mult_reduced)
