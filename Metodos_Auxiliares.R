#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
# Metodos auxiliares para analisis multicapa ------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------



# Matriz Laplaciana escalada ----------------------------------------------

calculate_scaled_laplacian <- function(mat) {

  nodes_degree <- rowSums(mat) + colSums(mat)

  mat_degree <- diag(x = nodes_degree, names = T,
                     nrow = length(nodes_degree),
                     ncol = length(nodes_degree))

  colnames(mat_degree) <- names(nodes_degree)

  c_norm <- 1/(2*sum(mat))

  # Calculamos Laplaciana escalada
  laplacian <- c_norm * (mat_degree - mat)

  return(laplacian)
}




# Entropia de Von Neumann -------------------------------------------------

calculate_von_neumann_entropy <- function(mat_lap) {
  eigenvalues_aux <- eigen(mat_lap)$values
  # Nos quedamos solo con los reales
  eigenvalues_aux <- eigenvalues_aux[Im(eigenvalues_aux) == 0]
  eigenvalues_aux <- Re(eigenvalues_aux)
  eigenvalues_aux <- eigenvalues_aux[eigenvalues_aux != 0]
  h_a <- (-1)*sum(eigenvalues_aux*log2(eigenvalues_aux))
  return(h_a)
}



# Divergencia de Jensen-Shannon -------------------------------------------
jensen_shannon_divergence <- function(mat_lap1, mat_lap2) {
  mixture <- (1/2)*(mat_lap1 + mat_lap2)
  entropy_mixture <- calculate_von_neumann_entropy(mixture)
  entropy_1 <- calculate_von_neumann_entropy(mat_lap1)
  entropy_2 <- calculate_von_neumann_entropy(mat_lap2)
  return(entropy_mixture - (1/2)*(entropy_1 + entropy_2))
}


# Greedy algorithm for multiplex reduction
greedy_reduction <- function(multilayer) {

    # Calculate laplacians
    laplacians <- map(.x = multilayer, .f = calculate_scaled_laplacian)
    # Extract the names of the layers
    layers <- setdiff(names(laplacians), "Agregada")

    # Create an empty matrix to store the Jensen-Shannon distances between
    # all pairs of layers
    mat_jsd <- matrix(0, nrow = length(layers), ncol = length(layers))

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
    min_jsd <- min(mat_jsd[mat_jsd > 0])

    # Get indices inside matrix to know which layers are
    indices <- which(mat_jsd == min_jsd, arr.ind = T)

    # Since the matriz might be symmetric, the previous command might return
    # a matrix, we keep only the rownames
    layers_to_merge <- rownames(indices) %>% sort() %>% unique()

    # We create the new multilayer system
    aux_new_multilayer <- multilayer[[layers_to_merge[1]]] + multilayer[[layers_to_merge[2]]]

    # Create the new multilayer system
    new_multilayer <- multilayer[setdiff(names(multilayer), layers_to_merge)]
    new_multilayer[[paste0(layers_to_merge[1], "+", layers_to_merge[2])]] <- aux_new_multilayer

    # Lets calculate new laplacians
    new_laplacians <- map(.x = new_multilayer, .f = calculate_scaled_laplacian)

    # New entropies
    new_entropies <- map(.x = new_laplacians, .f = calculate_von_neumann_entropy)

    # Entropy of new multilayer
    new_multilayer_entropy <- sum(unlist(new_entropies)) - new_entropies[["Agregada"]]

    # Create a new object that will be the result of this function
    result <- list(New_Multilayer = new_multilayer,
                   New_Total_Entropy = new_multilayer_entropy,
                   Aggregation = paste0(layers_to_merge[1], "+", layers_to_merge[2]))

    return(result)
}


# Recursive greedy reduction ----------------------------------------------

recursive_greedy_reduction <- function(multilayer) {

    # Primera reduccion
    multilayer_reduced_initial <- greedy_reduction(multilayer = multilayer)
    multilayer_reduced <- multilayer_reduced_initial

    # Define a condition to stop the iterations
    while(length(names(multilayer_reduced$New_Multilayer)) > 2) {
        aux <- greedy_reduction(multilayer = multilayer_reduced$New_Multilayer)
    }
}


# Agregacion optima de capas ----------------------------------------------

optimal_reduction <- function(multilayer) {

    if (!("Agregada" %in% names(multilayer))) stop("The aggregated layer of your multiplex should be called 'Agregada' in the list")

    # Lets find scaled laplacian matrices
    laplacians <- map(.x = multilayer, .f = calculate_scaled_laplacian)

    # With the laplacian matrices, we calculate the Von Neumann entropy
    entropies <- map(.x = laplacians, .f = calculate_von_neumann_entropy)

    # Let's calculate the entropy of the multiplex system
    multilayer_initial_entropy <- sum(unlist(entropies)) - entropies[["Agregada"]]

    # Since we will use it later, we store in a separate object the entropy
    # of the aggregated layer
    aggregated_layer_entropy <- entropies[["Agregada"]]

    # Let's create a vector that will keep the entropies of each step of
    # the aggregation process.
    entropies_by_aggregation <- vector(mode = "numeric")

    # Calculate initial entropy of the total disaggregated system
    entropies_by_aggregation <- multilayer_initial_entropy / (length(multilayer) - 1)

    # Add a name to distinguish the entropy of each aggregation step, to be able
    # to select the best one
    names(entropies_by_aggregation) <- "AllLayers"

    # Calculate all possible steps of the greedy aggregation algorithm





}
