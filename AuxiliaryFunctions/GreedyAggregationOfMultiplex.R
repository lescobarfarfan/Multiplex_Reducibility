#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
# Greedy algorithm to obtain optimal agregation of a multiplex ------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------

#' Obtain the optimal agregation of a multiplex network using the DeDomenico et. al (2015) method
#' @param multilayer A named list of matrices representing a multiplex network.
#' The aggregated matrix of the system is necessarily called 'Agregada'
#' @param character String to define the name of the aggregated layer
#' @param directed Logical value indicating if the layers are directed. Default FALSE.
#' @param weighted Logical value indicating if the layers are weighted. Default FALSE.
#' @return A list containing two items: the optimal aggregated network, and the
#' reducibility of the network
#' @export
greedy_reduction <- function(multilayer, name_agg, directed = F, weighted = F) {

    # Extract the names of the layers
    initial_layers <- setdiff(names(multilayer), name_agg)

    # Create initial values for a while loop that will perform the aggregation
    # of layers until optimal value is achieved.

    # Calculate laplacians and entropies per layer
    initial_laplacians <- map(.x = multilayer, .f = calculate_scaled_laplacian, directed, weighted)
    initial_entropies <- map(.x = initial_laplacians, .f = calculate_von_neumann_entropy)

    # Entropy of the aggregated system
    entropy_agg_sys <- initial_entropies[[name_agg]]

    # Entropy of the completeley differentiated system (all layers)
    initial_entropy_all_layers <- sum(unlist(initial_entropies[initial_layers]))
    initial_no_layers <- length(multilayer) - 1

    # Initial relative entropy
    initial_rel_entropy <- 1 - (initial_entropy_all_layers/initial_no_layers)*(1/entropy_agg_sys)

    # Initialize all objects that will store important values of each
    # aggregation obtained with the greedy algorithm
    aggregation_steps <- list(AllLayers = multilayer)
    relative_entropies_per_step <- vector(mode = "numeric")
    relative_entropies_per_step <- c(relative_entropies_per_step, initial_rel_entropy)
    names(relative_entropies_per_step) <- "AllLayers"

    # We initialize a while loop that will stop once the optmial aggregated
    # system is obtained
    rel_ent_control <- initial_rel_entropy
    agg_step <- 1

    # While loop
    while (rel_ent_control > 0 & length(aggregation_steps) < initial_no_layers) {
        if (agg_step == 1) {
            # Create an empty matrix to store the Jensen-Shannon distances between
            # all pairs of layers
            mat_jsd <- matrix(0, nrow = length(initial_layers),
                              ncol = length(initial_layers))

            # Assign names of the layers to identify them later for subsetting
            colnames(mat_jsd) <- initial_layers
            rownames(mat_jsd) <- initial_layers

            # Lets loop for all pairs of layers to calculate all J-S distances
            for (i in initial_layers) {
                for (j in initial_layers) {
                    lay1 <- initial_laplacians[[i]]
                    lay2 <- initial_laplacians[[j]]

                    aux_jsd <- jensen_shannon_divergence(mat_lap1 = lay1,
                                                         mat_lap2 = lay2)
                    mat_jsd[i, j] <- aux_jsd
                }
            }

            # lets find the smallest value (to merge such pair of matrices)
            # min_jsd <- abs(min(mat_jsd[mat_jsd != 0]))
            min_jsd <- min(mat_jsd[mat_jsd != 0])

            if (is.infinite(min_jsd)) {
                rel_ent_control <- 0
                break()
            }

            # Get indices inside matrix to know which layers are
            # indices <- which(abs(mat_jsd) == abs(min_jsd), arr.ind = T)
            indices <- which(mat_jsd == min_jsd, arr.ind = T)

            # Since the matriz might be symmetric, the previous command might return
            # a matrix, we keep only the rownames
            layers_to_merge <- rownames(indices) %>% sort() %>% unique()

            # We create the new multilayer system
            aux_new_multilayer <- multilayer[[layers_to_merge[1]]] + multilayer[[layers_to_merge[2]]]
            if (isFALSE(directed) & isFALSE(weighted)) {
                aux_new_multilayer[aux_new_multilayer > 1] <- 1
            } else if (isFALSE(weighted)) {
                aux_new_multilayer[aux_new_multilayer > 1] <- 1
            }

            # Create the new multilayer system
            new_multilayer <- multilayer[setdiff(names(multilayer), layers_to_merge)]
            new_multilayer[[paste0(layers_to_merge[1], "+", layers_to_merge[2])]] <- aux_new_multilayer

            # Aggregate to the list keeping track of all aggregations steps
            aggregation_steps[[paste0(layers_to_merge[1], "+", layers_to_merge[2])]] <- new_multilayer

            # Calculations of all new metrics for new aggregated system
            layers <- setdiff(names(new_multilayer), name_agg)
            laplacians <- map(.x = new_multilayer, .f = calculate_scaled_laplacian)
            entropies <- map(.x = laplacians, .f = calculate_von_neumann_entropy)

            entropy_all_layers <- sum(unlist(entropies[layers]))
            no_layers <- length(layers)
            rel_entropy <- 1 - (entropy_all_layers/no_layers)*(1/entropy_agg_sys)

            relative_entropies_per_step <- c(relative_entropies_per_step,
                                             rel_entropy)
            names(relative_entropies_per_step)[agg_step + 1] <- paste0(layers_to_merge[1],
                                                                       "+",
                                                                       layers_to_merge[2])
            agg_step <- 2
            rel_ent_control <- rel_entropy

        } else {

            # Continuamos la agregacion de capas
            # Create an empty matrix to store the Jensen-Shannon distances between
            # all pairs of layers
            mat_jsd <- matrix(0, nrow = length(layers),
                              ncol = length(layers))

            # Assign names of the layers to identify them later for subsetting
            colnames(mat_jsd) <- layers
            rownames(mat_jsd) <- layers

            # Lets loop for all pairs of layers to calculate all J-S distances
            for (i in layers) {
                for (j in layers) {
                    lay1 <- laplacians[[i]]
                    lay2 <- laplacians[[j]]

                    aux_jsd <- jensen_shannon_divergence(mat_lap1 = lay1,
                                                         mat_lap2 = lay2)
                    mat_jsd[i, j] <- aux_jsd
                }
            }

            # lets find the smallest value (to merge such pair of matrices)
            # min_jsd <- abs(min(mat_jsd[mat_jsd != 0]))
            min_jsd <- min(mat_jsd[mat_jsd != 0])

            if (is.infinite(min_jsd)) {
                rel_ent_control <- 0
                break()
            }

            # Get indices inside matrix to know which layers are
            # indices <- which(abs(mat_jsd) == abs(min_jsd), arr.ind = T)
            indices <- which(mat_jsd == min_jsd, arr.ind = T)

            # Since the matriz might be symmetric, the previous command might return
            # a matrix, we keep only the rownames
            layers_to_merge <- rownames(indices) %>% sort() %>% unique()

            # We create the new multilayer system
            aux_new_multilayer <- new_multilayer[[layers_to_merge[1]]] + new_multilayer[[layers_to_merge[2]]]
            if (isFALSE(directed) & isFALSE(weighted)) {
                aux_new_multilayer[aux_new_multilayer > 1] <- 1
            } else if (isFALSE(weighted)) {
                aux_new_multilayer[aux_new_multilayer > 1] <- 1
            }

            # Create the new multilayer system
            new_multilayer <- new_multilayer[setdiff(names(new_multilayer),
                                                     layers_to_merge)]
            new_multilayer[[paste0(layers_to_merge[1], "+", layers_to_merge[2])]] <- aux_new_multilayer

            # Aggregate to the list keeping track of all aggregations steps
            aggregation_steps[[paste0(layers_to_merge[1], "+", layers_to_merge[2])]] <- new_multilayer

            # Calculations of all new metrics for new aggregated system
            layers <- setdiff(names(new_multilayer), name_agg)
            laplacians <- map(.x = new_multilayer, .f = calculate_scaled_laplacian)
            entropies <- map(.x = laplacians, .f = calculate_von_neumann_entropy)

            entropy_all_layers <- sum(unlist(entropies[layers]))
            no_layers <- length(layers)
            rel_entropy <- 1 - (entropy_all_layers/no_layers)*(1/entropy_agg_sys)

            relative_entropies_per_step <- c(relative_entropies_per_step,
                                             rel_entropy)
            names(relative_entropies_per_step)[agg_step + 1] <- paste0(layers_to_merge[1],
                                                                       "+",
                                                                       layers_to_merge[2])
            agg_step <- agg_step + 1
            rel_ent_control <- rel_entropy

        }
    }

    # Once the algorithm is finish, we calculate the reducibility and return
    # the optimal agregation
    # max_rel_agg <- which(relative_entropies_per_step == max(relative_entropies_per_step[-"AllLayers")]))
    max_rel_agg <- sort(relative_entropies_per_step, decreasing = T)[2]
    name_optimal_agregation <- names(max_rel_agg)


    optimal_agregation <- aggregation_steps[[name_optimal_agregation]]

    reducibility <- ((length(multilayer) - 1) - (length(optimal_agregation) - 1)) / (length(multilayer) - 1)

    output <- list(OptimalAgg = optimal_agregation,
                   Reducibility = reducibility)

    return(output)
}