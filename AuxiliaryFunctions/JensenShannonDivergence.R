#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
# Jensen-Shannon divergence between two graphs ----------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------


#' Obtain the divergence betwenn two laplacian (density) matrices
#'
#' @param mat_lap1 A square matrix, the laplacian of a graph
#' @param mat_lap2 A square matrix, the laplacian of a graph
#' @return The divergence between two graphs
#' @export
jensen_shannon_divergence <- function(mat_lap1, mat_lap2) {
    mixture <- (1/2)*(mat_lap1 + mat_lap2)
    entropy_mixture <- calculate_von_neumann_entropy(mixture)
    entropy_1 <- calculate_von_neumann_entropy(mat_lap1)
    entropy_2 <- calculate_von_neumann_entropy(mat_lap2)
    return(entropy_mixture - (1/2)*(entropy_1 + entropy_2))
}