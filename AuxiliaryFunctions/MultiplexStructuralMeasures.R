#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
# Multiplex structural measures -------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------
#   -----------------------------------------------------------------------


#' Function to calculate the multiplex degree of nodes
#' 
#' @param multilayer A multilayer network in the form of a list
#' @param directed Boolean value to indicate if the networks is directed
#' @return A data frame containing the degree of each node in each layer
#' @importFrom purrr map
#' @export
calculate_multiplex_degree <- function(multilayer, directed = FALSE) {
  
  if (directed) {
    degrees <- map(.x = multilayer, .f = function(mat) return(rowSums(mat) + colSums(mat)))
  } else {
    degrees <- map(.x = multilayer, .f = colSums)
  }
  degrees_df <- bind_cols(degrees)
  degrees_df$Node <- names(degrees[[1]])
  degrees_df <- degrees_df %>% select(Node, 1:(ncol(degrees_df)-1))
  degrees_df <- degrees_df %>% mutate(OverlappingDeg = rowSums(degrees_df[2:ncol(degrees_df)]) - Agregada)
  return(degrees_df)
}


#' Function to calculate the participation 