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
#' @return A data frame containing the degree of each node in each layer. If
#' the network is directed, it will return a list containing three elements:
#' the in degree, the out degree, and the total degree. If the system is weighted
#' then the strenght will be calculated.
#' @importFrom purrr map
#' @export
calculate_multiplex_degree <- function(multilayer, directed = FALSE,
                                       weighted = FALSE) {

  if (isTRUE(directed) & isTRUE(weighted)) {
    degree_out <- map(.x = multilayer, .f = function(mat) return(rowSums(mat)))
    degree_in <- map(.x = multilayer, .f = function(mat) return(colSums(mat)))
    degree_tot <- map(.x = multilayer, .f = function(mat) return(rowSums(mat) + colSums(mat)))

    deg_in_df <- bind_cols(degree_in)
    deg_in_df$Node <- names(degree_in[[1]])
    deg_in_df <- deg_in_df %>%
      select(Node, 1:(ncol(deg_in_df) - 1))
    deg_in_df <- deg_in_df %>%
      mutate(OverlappingDeg_In = rowSums(deg_in_df[2:ncol(deg_in_df)]) - Agregada)

    degree_out_df <- bind_cols(degree_out)
    degree_out_df$Node <- names(degree_in[[1]])
    degree_out_df <- degree_out_df %>%
      select(Node, 1:(ncol(degree_out_df) - 1))
    degree_out_df <- degree_out_df %>%
      mutate(OverlappingDeg_Out = rowSums(degree_out_df[2:ncol(degree_out_df)]) - Agregada)

    degree_tot_df <- bind_cols(degree_tot)
    degree_tot_df$Node <- names(degree_tot[[1]])
    degree_tot_df <- degree_tot_df %>%
      select(Node, 1:(ncol(degree_tot_df) - 1))
    degree_tot_df <- degree_tot_df %>%
      mutate(OverlappingDeg_Tot = rowSums(degree_tot_df[2:ncol(degree_tot_df)]) - Agregada)

    degrees_df <- list(Deg_In = deg_in_df, Deg_Out = degree_out_df, Deg_Tot = degree_tot_df)
    return(degrees_df)

  } else if (isTRUE(directed) & isFALSE(weighted)) {
    degree_out <- map(.x = multilayer, .f = function(mat) {
      mat[mat > 0] <- 1
      return(rowSums(mat))
    })
    degree_in <- map(.x = multilayer, .f = function(mat) {
      mat[mat > 0] <- 1
      return(colSums(mat))
    })
    degree_tot <- map(.x = multilayer, .f = function(mat) {
      mat[mat > 0] <- 1
      return(rowSums(mat) + colSums(mat))
    })

    deg_in_df <- bind_cols(degree_in)
    deg_in_df$Node <- names(degree_in[[1]])
    deg_in_df <- deg_in_df %>%
      select(Node, 1:(ncol(deg_in_df) - 1))
    deg_in_df <- deg_in_df %>%
      mutate(OverlappingDeg_In = rowSums(deg_in_df[2:ncol(deg_in_df)]) - Agregada)

    degree_out_df <- bind_cols(degree_out)
    degree_out_df$Node <- names(degree_in[[1]])
    degree_out_df <- degree_out_df %>%
      select(Node, 1:(ncol(degree_out_df) - 1))
    degree_out_df <- degree_out_df %>%
      mutate(OverlappingDeg_Out = rowSums(degree_out_df[2:ncol(degree_out_df)]) - Agregada)

    degree_tot_df <- bind_cols(degree_tot)
    degree_tot_df$Node <- names(degree_tot[[1]])
    degree_tot_df <- degree_tot_df %>%
      select(Node, 1:(ncol(degree_tot_df) - 1))
    degree_tot_df <- degree_tot_df %>%
      mutate(OverlappingDeg_Tot = rowSums(degree_tot_df[2:ncol(degree_tot_df)]) - Agregada)

    degrees_df <- list(Deg_In = deg_in_df, Deg_Out = degree_out_df, Deg_Tot = degree_tot_df)
    return(degrees_df)

  } else if (isFALSE(directed) & isTRUE(weighted)) {
    degrees <- map(.x = multilayer, .f = colSums)
    degrees_df <- bind_cols(degrees)
    degrees_df$Node <- names(degrees[[1]])
    degrees_df <- degrees_df %>%
      select(Node, 1:(ncol(degrees_df) - 1))
    degrees_df <- degrees_df %>%
      mutate(OverlappingDeg = rowSums(degrees_df[2:ncol(degrees_df)]) - Agregada)
    return(degrees_df)
  } else if (isFALSE(directed) & isTRUE(weighted)) {
    degrees <- map(.x = multilayer, .f = colSums)
    degrees_df <- bind_cols(degrees)
    degrees_df$Node <- names(degrees[[1]])
    degrees_df <- degrees_df %>%
      select(Node, 1:(ncol(degrees_df) - 1))
    degrees_df <- degrees_df %>%
      mutate(OverlappingDeg = rowSums(degrees_df[2:ncol(degrees_df)]) - Agregada)
    return(degrees_df)

  }

}
