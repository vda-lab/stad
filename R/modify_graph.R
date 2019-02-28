
# ----------------------
# Table of contents
# ----------------------
# modify_graph
# ----------------------

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
#' Modify graph
#'
#' Update a STAD graph with the a new number of edges.
#'
#' @param x \code{stad} object
#' @param number_edges string variable name for the first dimension of the lens. Default 'from_x'.
#'
#' @export
#' @return Returns a \code{stad} class (or \code{list}) with the following items:
#' \itemize{
#' \item graph
#' \item number_edges
#' \item trace
#' \item ordered_distances
#' }
modify_graph <- function(x, number_edges){

  if (class(x) != "stad") stop("It is not a stad object.")

  # Genearte the array of vertices
  vertices_names <- as.numeric(igraph::vertex_attr(x$graph)$name)

  # Return process
  list_to_return <- list(
    graph = generate_graph(ordered_distances = x$ordered_distances,
                           number_edges = number_edges,
                           vertices_names = vertices_names
                           ),
    number_edges = number_edges,
    trace = NA,
    ordered_distances = x$ordered_distances
  )
  class(list_to_return) <- "stad"

  return(list_to_return)
}
