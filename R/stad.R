
# ----------------------
# Table of contents
# ----------------------
# stad
# ----------------------

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
#' STAD algorithm
#'
#' The function computes STAD with a given distance matrix and/or filter values.
#'
#' @param distance_matrix dist object containing a distance matrix.
#' @param filter_values numeric array or two-dimensional data.frame with the lenses to use in STAD.
#' @param num_intervals number or array defining the number of intervals to consider in the filters.
#' @param metric string or array defining the metrics supported ("polar" or "euclidean").
#' @param penalty numeric penalty value used during the optimization loop. It limits the number of iterations. Default = 0.
#' @param random_factor integer factor that controls the variability of the next iteration on the optimization algorithm. Default = 10000. Higher values generates more diverse exploration.
#' @param iterations_inner_loop integer number of evaluations for each iteration. The number of iterations allows knowing the close values around but increase the total number of iterations.
#' @param two_mst boolean indicating if the MST is build using the two-step MST. It only applies when \code{filter_values} are defined.
#'
#' @export
#' @return Returns a \code{stad} class (or \code{list}) with the following items:
#' \itemize{
#' \item graph
#' \item number_edges
#' \item trace
#' \item ordered_distances
#' }
#' @examples
#' # Iris dataset
#' data(iris)
#' iris_distance <- dist(iris[,1:4])
#' iris_stad <- stad(iris_distance)
#' plot_graph(iris_stad)
#' plot_trace(iris_stad)
#' # Circles dataset
#' data(circles)
#'
#' library(magrittr)
#' library(ggplot2)
#'
#' circles %>%
#'   ggplot(aes(x,y, color = lens)) +
#'   geom_point()
#'
#' circles_distance <- dist(circles[,c("x", "y")])
#'
#' ## STAD without lens
#' set.seed(10)
#' circles_nolens <- stad(circles_distance)
#' plot_graph(circles_nolens, layout = igraph::layout_with_kk )
#'
#' ## STAD with lens
#' set.seed(10)
#' circles_lens <- stad(circles_distance, filter_values = circles$lens, num_intervals = 5)
#' plot_graph(circles_lens, layout = igraph::layout_with_kk )
#'
stad <- function( distance_matrix,
                  filter_values = NULL,
                  num_intervals = NULL,
                  metric = NULL,
                  penalty = 0,
                  random_factor = 10000,
                  iterations_inner_loop = 10,
                  two_mst = FALSE){

  if (is.null(filter_values)) {

    stad_without_lens( distance_matrix,
                       penalty = penalty,
                       random_factor = random_factor,
                       iterations_inner_loop = iterations_inner_loop)

  } else if(!is.null(filter_values)) {

    stad_with_lens ( distance_matrix,
                     filter_values = filter_values,
                     num_intervals = num_intervals,
                     metric = metric,
                     penalty = penalty,
                     random_factor = random_factor,
                     iterations_inner_loop = iterations_inner_loop,
                     two_mst = two_mst)

  }


}
