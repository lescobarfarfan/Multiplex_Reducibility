#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
# "main" para calcular metricas de reducibilidad --------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------



# Paquetes ----------------------------------------------------------------

library(tidyverse)
library(Matrix)



# corremos codigos previos al analisis de reducibilidad -------------------
source("C:/Users/K15523/Desktop/Multiplex_Reducibility/Metodos_Auxiliares.R")


# Cargamos el objeto multicapa (lista de listas) --------------------------

multicapa_temporal <- readRDS("C:/Users/K15523/Desktop/Multiplex_Reducibility_Auxiliares/Resultados/MulticapaTemporal_Ad.Rds")



# Lineas temporales para una sola fecha -----------------------------------

multicapa_test <- multicapa_temporal[[449]]





# # Calculate Laplacian matrices
# laplacianas <- map(.x = multicapa_test, calcula_laplaciana_escalada)
#
# entropias <- map(.x = laplacianas, calcula_entropia_von_neumann)
#
# entropia_multicapa <- sum(unlist(entropias)) - entropias[["Agregada"]]
#
# entropia_agregada <- entropias[["Agregada"]]
#
# entropias_multicapa <- vector(mode = "numeric")
# entropias_multicapa <- entropia_multicapa/(length(laplacianas)-1)
# names(entropias_multicapa) <- "AllLayers"
#
# # test --------------------------------------------------------------------
#
#
#
# # Pair-wise Jensen-Shannon divergence -------------------------------------
# capas <- c("CVT", "Reportos", "Unsecured")
# matriz_jensen_shannon_pairs <- matrix(data = 0, nrow = length(capas),
#                                       ncol = length(capas))
# colnames(matriz_jensen_shannon_pairs) <- capas
# rownames(matriz_jensen_shannon_pairs) <- capas
#
#
# # Para cada par de capas, calculamos la divergencia de Jensen-Shannon
# for (i in capas) {
#   for (j in capas) {
#     mat_lap1 <- laplacianas[[i]]
#     mat_lap2 <- laplacianas[[j]]
#     aux_jsd <- jensen_shannon_div(mat_lap1 = mat_lap1, mat_lap2 = mat_lap2)
#     matriz_jensen_shannon_pairs[i, j] <- aux_jsd
#   }
# }
#
# # Con la matriz, hacemos el siguiente paso del algoritmo: sumar las matrices
# # con el valor mas pequenio de la divergencia Jensen-Shannon
# min_jsd <- min(matriz_jensen_shannon_pairs[matriz_jensen_shannon_pairs > 0])
# indices <- which(matriz_jensen_shannon_pairs == min_jsd, arr.ind = T)
# layers_to_merge <- rownames(indices) %>% unique() %>% sort()
#
# aux_new_multilayer <- multicapa_test[[layers_to_merge[1]]] +
#                       multicapa_test[[layers_to_merge[[2]]]]
#
# new_multilayer <- multicapa_test[setdiff(names(multicapa_test), layers_to_merge)]
#
#
# new_multilayer[[paste0(layers_to_merge[1], "+", layers_to_merge[2])]] <- aux_new_multilayer
#
# laplacianas_new <- map(.x = new_multilayer, .f = calcula_laplaciana_escalada)
# entropias_new <- map(.x = laplacianas_new, .f = calcula_entropia_von_neumann)
#
#
# entropia_multicapa_new <- sum(unlist(entropias_new)) - entropias_new[["Agregada"]]
#
#
# entropias_multicapa <- c(entropias_multicapa, entropia_multicapa_new/(length(laplacianas_new)-1))
# names(entropias_multicapa)[2] <- paste0(layers_to_merge[1], "+", layers_to_merge[2])
#
# distinguishabilities <- 1 - entropias_multicapa/entropia_agregada
# names(distinguishabilities) <- names(entropias_multicapa)
# which(distinguishabilities == max(distinguishabilities))
