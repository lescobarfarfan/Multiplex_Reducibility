#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
# Script to calculate Jaccard Index between two layers --------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------

#' Function to calculate the similarity between two layers with the same set
#' of nodes using the Jaccard index
#' @param matrix1 An adjacency matrix representing a network.
#' @param matrix2 An adjacency matrix representing a network. 
#' @return Returns a numeric value representing the level of similarity between
#' two networks.
#' @export
similarity_between_layers <- function(matrix1, matrix2) {
  if ((nrow(matrix1) != nrow(matrix2)) | (ncol(matrix1) != ncol(matrix2))) {
    stop("Matrices must be of the same size")
  }
  if (any(matrix1 > 1) | any(matrix2 > 1)) {
    stop("Only adjacency matrix are allowed")
  }
  mat_aux <- matrix1 + matrix2
  mat_aux[mat_aux > 1] <- 1
  similarity <- sum(mat_aux)/(sum(matrix1) + sum(matrix2))
  return(similarity)
}