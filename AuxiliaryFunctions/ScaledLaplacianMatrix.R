#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
# Function to calculate the scaled laplacian matrix of a graph ------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------


#' Calculate the Laplacian function of a directed and weighted graph
#'
#' @param mat A square adjacency matrix of a directed (and possible weighted) graph
#' @return A squared matrix representing the Laplacian of a directed (and
#' possibly weighted) graph.
laplacian_directed <- function(mat, weighted = F) {

    # First, we have to turn the matrix into a data frame
    mat_df <- mat %>%
        as.data.frame() %>%
        tbl_df()
    mat_df$Inst <- colnames(mat)
    mat_df <- mat_df %>%
        select(Inst, 1:(ncol(mat_df) - 1))
    mat_df <- mat_df %>%
        gather(key = Cont, value = Exp, -Inst)

    # Delete self loops and exposures equalt to zero
    mat_df <- mat_df %>%
        filter(Exp > 0,
               Inst != Cont)

    #Ad an ID to each edge
    mat_df$ID_E <- seq(1, nrow(mat_df), 1)

    # Extrac auxiliary vectors
    from <- mat_df$Inst
    to <- mat_df$Cont
    edges <- mat_df$ID_E

    rownames(mat) <- colnames(mat)

    # Extract values to create dimension of insidence matrix
    n <- length(colnames(mat))
    m <- nrow(mat_df)

    # Create empty matrix to store incidence matrix values
    incidence_mat <- matrix(data = 0, nrow = n, ncol = m)
    rownames(incidence_mat) <- rownames(mat)
    colnames(incidence_mat) <- mat_df$ID_E

    # Fill the incidence matrix
    for (i in 1:nrow(mat_df)) {
        from_aux <- from[i]
        to_aux <- to[i]
        edge_aux <- as.character(edges[i])
        row_from_indx_aux <- which(rownames(incidence_mat) == from_aux)
        row_to_indx_aux <- which(rownames(incidence_mat) == to_aux)
        col_indx_aux <- which(colnames(incidence_mat) == edge_aux)
        incidence_mat[row_from_indx_aux, col_indx_aux] <- 1
        incidence_mat[row_to_indx_aux, col_indx_aux] <- -1
    }

    # Calculate laplacian
    if (weighted) {
        # Create a diagonal matrix with the weights
        exposures <- mat_df$Exp
        exposures_mat <- diag(x = exposures, nrow = m, ncol = m)

        laplacian <- incidence_mat %*% exposures_mat %*% t(incidence_mat)
    } else {
        laplacian <- incidence_mat %*% t(incidence_mat)
    }

    return(laplacian)

}



#' Calculate the scaled laplacian function of a graph
#'
#' @param mat A square adjacency matrix of a graph
#' @param directed Boolean value to determine if mat is directed, in which case
#' @param weighted Boolean value to determine if mat is weighted, in which case
#' the total weighted degree is calculate
#' the total degree is calculated
#' @return A square matrix representing the scales laplacian matrix of a graph
#' @export
calculate_scaled_laplacian <- function(mat, directed = F, weighted = F) {

    # if (isTRUE(directed) & isTRUE(weighted)) {
    #     nodes_degree <- rowSums(mat) + colSums(mat)
    # } else if (isTRUE(directed) & isFALSE(weighted)) {
    #     nodes_degree <- rowSums(mat) + colSums(mat)
    # } else if (isFALSE(directed) & isTRUE(weighted)) {
    #     nodes_degree = rowSums(mat)
    # } else {
    #     nodes_degree = rowSums(mat)
    # }

    if (isTRUE(directed) & isTRUE(weighted)) {
        laplacian <- laplacian_directed(mat = mat, weighted = T)
        # Scaling
        c_norm <- 1/(2*sum(mat))

        # Calculamos Laplaciana escalada
        laplacian <- c_norm * laplacian

    } else if (isTRUE(directed) & isFALSE(weighted)) {
        laplacian <- laplacian_directed(mat = mat, weighted = F)
        # Scaling
        c_norm <- 1/(2*sum(mat))

        # Calculamos Laplaciana escalada
        laplacian <- c_norm * laplacian

    } else if (isFALSE(directed) & isTRUE(weighted)) {
        nodes_degree <- rowSums(mat)
        mat_degree <- diag(x = nodes_degree, names = T,
                           nrow = length(nodes_degree),
                           ncol = length(nodes_degree))

        colnames(mat_degree) <- names(nodes_degree)

        c_norm <- 1/(2*sum(mat))

        # Calculamos Laplaciana escalada
        laplacian <- c_norm * (mat_degree - mat)
    } else if (isFALSE(directed) & isFALSE(weighted)) {
        nodes_degree <- rowSums(mat)
        mat_degree <- diag(x = nodes_degree, names = T,
                           nrow = length(nodes_degree),
                           ncol = length(nodes_degree))

        colnames(mat_degree) <- names(nodes_degree)

        c_norm <- 1/(2*sum(mat))

        # Calculamos Laplaciana escalada
        laplacian <- c_norm * (mat_degree - mat)
    }

    return(laplacian)
}
