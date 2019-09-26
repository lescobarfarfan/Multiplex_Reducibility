#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
# Function to calculate the scaled laplacian matrix of a graph ------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------


#' Calculate the scaled laplacian function of a graph
#'
#' @param mat A square adjacency matrix of a graph
#' @return A square matrix representing the scales laplacian matrix of a graph
#' @export
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