#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
# Script to calculate Von-Neumann entropy ---------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------


#' Calculate von Neumann entropy of a graph
#'
#' @param mat_lap A square matrix representing the laplacian of a graphh
#' @param tol A double representing the error tolerance for floating point error
#' @return A double representing the Von Neumann entropy of a graph
#' @export
calculate_von_neumann_entropy <- function(mat_lap, tol = 10^-6) {
    eigenvalues_aux <- eigen(mat_lap)$values
    eigenvalues_aux <- eigenvalues_aux[Im(eigenvalues_aux) == 0]
    eigenvalues_aux <- Re(eigenvalues_aux)
    eigenvalues_aux <- eigenvalues_aux[abs(eigenvalues_aux - 0) > tol]
    h_a <- (-1)*sum(eigenvalues_aux*log2(eigenvalues_aux))
    return(h_a)
}


