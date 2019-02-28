
# ----------------------
# Table of contents
# ----------------------
# plot_trace
# plot_graph
# ----------------------

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
#' Plot trace
#'
#' Line-chart of the intermediate evaluations during the optimizations process.
#'
#' @param x \code{stad} object.
#' @importFrom rlang .data
#' @export
#' @examples
#' # Iris dataset
#' data(iris)
#' iris_distance <- dist(iris[,1:4])
#' iris_stad <- stad(iris_distance)
#' plot_graph(iris_stad)
#' plot_trace(iris_stad)
#'
plot_trace <- function(x){

  if( class(x) != "stad") stop("It is not a list object.")

  # Define x range
  range_x <- range(x$trace$iteration)

  # Adapt labels for Max correlation, Selected number of edges or Other evaluations
  x$trace <- x$trace %>%
    dplyr::mutate(
      type = ifelse(
        .data$correlation == max( x$trace$correlation ),
        "Max. correlation",
        "Other evaluations"
      ),
      type = ifelse(
        .data$type == "Other evaluations" & .data$iteration == x$number_edges,
        "Selected number of edges",
        .data$type
      ),
      type = factor(.data$type,  levels = c("Other evaluations", "Max. correlation", "Selected number of edges"))
    )

  x$trace %>%
    ggplot2::ggplot(ggplot2::aes(.data$iteration, .data$correlation, color  = .data$type)) +
    ggplot2::geom_line(colour = "#252525") +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(
      values = c("#252525", "#fc4e2a", "#fe9929")
    ) +
    ggplot2::geom_point(
      data = x$trace %>%
        dplyr::filter(.data$type == "Max. correlation"), fill = "#fc4e2a", shape = 21, color = "#252525", size = 2.5) +
    ggplot2::geom_point(
      data = x$trace %>%
        dplyr::filter(.data$type == "Selected number of edges"), fill = "#fe9929", shape = 21, color = "#252525", size = 2.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="bottom") +
    ggplot2::scale_x_continuous(breaks = round(seq(range_x[1], range_x[2], range_x[2] /10) ) ) +
    ggplot2::scale_y_continuous(limits=c(0, 1), breaks = seq(0,1, 0.1)) +
    ggplot2::labs(x = "Num. edges",
                  y = "Objective Function (Correlation)",
                  title = "Intermediate results during STAD process",
                  color = "",
                  caption = "Generated with plot_trace")
}


# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
#' Plot graph
#'
#' Graph plot for the resulting STAD graph.
#'
#' @param x \code{stad} object.
#' @param layout igraph layout. Default igraph::layout_nicely
#' @param vertex.color Vertex color for the plot. Default black.
#' @param vertex.frame.color Vertex color contour. Default grey.
#' @param vertex.size Vertex size. Default 3
#' @param edge.width Edge width. Default 0.5
#' @export
#' @examples
#' # Iris dataset
#' data(iris)
#' iris_distance <- dist(iris[,1:4])
#' iris_stad <- stad(iris_distance)
#' plot_graph(iris_stad)
#' plot_trace(iris_stad)
#'
plot_graph <- function(x,
                       layout = igraph::layout_nicely,
                       vertex.color = "black",
                       vertex.frame.color = "grey",
                       vertex.size = 3,
                       edge.width = 0.5){
  igraph::plot.igraph(x$graph,
       layout = layout,
       vertex.color = vertex.color,
       vertex.frame.color = vertex.frame.color,
       vertex.size = vertex.size,
       edge.width = edge.width,
       vertex.label=NA)
}
