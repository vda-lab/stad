
# ----------------------
# Table of contents
# ----------------------
# stad_without_lens
# stad_evaluation
# modify_distance_matrix
# generate_graph
# modify_output
# ----------------------

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
#' STAD without lens
#'
#' The function computes STAD without lens using a given distance matrix.
#'
#' @param distance_matrix dist object containing a distance matrix.
#' @param penalty numeric penalty value used during the optimization loop. It limits the number of iterations. Default = 0.
#' @param random_factor integer factor that controls the variability of the next iteration on the optimization algorithm. Default = 10000. Higher values generates more diverse exploration.
#' @param iterations_inner_loop integer number of evaluations for each iteration. The number of iterations allows knowing the close values around but increase the total number of iterations.
#'
#' @return Returns a \code{stad} class (or \code{list}) with the following items:
#' \itemize{
#' \item graph
#' \item number_edges
#' \item trace
#' \item ordered_distances
#' }
stad_without_lens <- function( distance_matrix,
                               penalty = 0,
                               random_factor = 10000,
                               iterations_inner_loop = 10,
                               ratio = FALSE){

  # Distance matrix transformation
  modified_distance_matrix <- modify_distance_matrix(distance_matrix)
  vertices_names <- modified_distance_matrix$vertices_names
  distance_matrix <- modified_distance_matrix$distance_matrix
  distance_matrix_df <- modified_distance_matrix$distance_matrix_df

  # Distance matrix is transformed into a graph (igraph)
  distance_matrix_graph <- igraph::graph_from_data_frame(distance_matrix_df,
                                                         directed = FALSE,
                                                         vertices = vertices_names)

  # STAD iteration
  list_to_return <- stad_evaluation(distance_matrix = distance_matrix,
                                    distance_matrix_df = distance_matrix_df,
                                    vertices_names = vertices_names,
                                    mst_graph = igraph::mst(distance_matrix_graph,
                                                            weights = distance_matrix_df$value),
                                    penalty = penalty,
                                    random_factor = random_factor,
                                    iterations_inner_loop = iterations_inner_loop,
                                    ratio = ratio)

  return(list_to_return)
}


# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
#' STAD evaluation
#'
#' Iterative process in STAD using Simulated Annealing algorithm
#'
#' @param distance_matrix Distance matrix as \code{dist} format.
#' @param distance_matrix_df Distance matrix as \code{data.frame} format.
#' @param vertices_names numeric array with names of nodes (\code{colnames(distance_matrix)})
#' @param mst_graph Minimum Spanning Tree graph from distance_matrix
#' @param penalty numeric penalty value used during the optimization loop. It limits the number of iterations. Default = 0.
#' @param random_factor integer factor that controls the variability of the next iteration on the optimization algorithm. Default = 10000. Higher values generates more diverse exploration.
#' @param iterations_inner_loop integer number of evaluations for each iteration. The number of iterations allows knowing the close values around but increase the total number of iterations.
#' @importFrom rlang .data
#'
#' @return Returns a \code{stad} class (or \code{list}) with the following items:
#' \itemize{
#' \item graph
#' \item number_edges
#' \item trace
#' \item ordered_distances
#' }
stad_evaluation <- function( distance_matrix = NULL,
                             distance_matrix_df = NULL,
                             vertices_names = NULL,
                             mst_graph = NULL,
                             penalty = NULL,
                             random_factor = NULL,
                             iterations_inner_loop = NULL,
                             ratio = FALSE){

  # Computing Minimum Spanning Tree from igraph package (MST)
  mst_matrix <- igraph::as_data_frame(mst_graph) %>%
    dplyr::mutate(from = as.numeric(.data$from),
                  to = as.numeric(.data$to),
                  mst = 1) %>%
    dplyr::select(.data$from, .data$to, .data$mst)

  # Merge MST information with Distance Matrix
  spanning_tree <- dplyr::full_join(
    distance_matrix_df, mst_matrix, by = c("from", "to")
  ) %>%
    dplyr::arrange(.data$value) %>%
    dplyr::mutate(mst = ifelse(is.na(.data$mst), 0, .data$mst)) %>%
    dplyr::group_by(.data$mst) %>%
    dplyr::mutate(order = dplyr::row_number()) %>%
    dplyr::ungroup()

  # Max number of edges in the network
  max_iterations <- nrow(spanning_tree %>% dplyr::filter(.data$mst == 0))

  # ---
  # Constant value used in the penalty
  constant_size <- sqrt( 1/max_iterations )

  # 0.001 * x * const = 1 - evaluation(0)
  # x = (1 - evaluation(0)) / (pena * const)
  # ----

  # Define distance matrix array
  # Adding constant values for the the distance matrix
  distance_matrix_array <- as.numeric(distance_matrix)

  # Evaluation function to be optimized
  evaluation <- function (i, verbose = TRUE){
    # Select the the number of edges to be evaluated
    iteration <- spanning_tree %>% dplyr::filter( order <= i | .data$mst == 1 )
    # Ratio similarity-distance
    ratio <- sum(1-iteration$value) / ( 1 + sum(iteration$value) )
    # Transformation into a graph
    iteration_graph <- igraph::graph_from_data_frame(iteration, directed = FALSE, vertices = vertices_names)
    # Computing shortest path matrix
    iteration_dist <- igraph::distances(iteration_graph)
    # Add unit-distance matrix to compute the correlation
    comparsion_matrix <- cbind(distance_matrix_array, as.numeric(iteration_dist))
    colnames(comparsion_matrix) <- c("a", "b")
    comparsion_matrix <- as.data.frame(comparsion_matrix) %>% dplyr::filter(!is.na(.data$a) & !is.na(.data$b))
    # Computing Pearson correlation
    # The filter remove NA's value from stad with lens
    correlation <- stats::cor( comparsion_matrix$a, comparsion_matrix$b)
    # Print intermediate result to be captured
    if (verbose) print( c( i, correlation) )
    #
    if (ratio == TRUE) {
      obj_fun <- ratio * correlation - ( penalty * i *  constant_size )
    } else {
      obj_fun <- correlation - ( penalty * i *  constant_size )
    }
    return( obj_fun )
  }

  # ----
  upper_max_iterations <- (1 - evaluation(0, verbose = FALSE)) / (penalty * constant_size)
  # ----

  # Function that determines the next iteration
  next_iteration <- function(para_0, fun_length, rf, temp = NA){
    # Integer within the random_factor (rf) range
    sample_value = sample.int(rf, fun_length, replace = TRUE)
    # Binomial function. It returns +1 or -1
    binom = ((stats::rbinom(fun_length, 1, 0.5) * -2) + 1)
    return (para_0 +  sample_value * binom)
  }

  # Execution SA and capture the intermediate results.
  intermediate_results <- utils::capture.output(

    optimization::optim_sa(evaluation, 0, maximization = TRUE, trace = FALSE,
                           lower = 0,
                           upper = ifelse(
                             penalty == 0 |  upper_max_iterations > max_iterations,
                             max_iterations,
                             upper_max_iterations
                           ),
                           control = list(vf = next_iteration,
                                          rf = random_factor,
                                          ac_acc = 0.005,
                                          maxgood = iterations_inner_loop,
                                          stopac = 1,
                                          nlimit = iterations_inner_loop)
    )

  )

  modified_output <- modify_output(intermediate_results)

  # Print num of iterations
  print(
    paste("Network with", modified_output$max_value,
          "edges of", max_iterations, "and", length(vertices_names), "vertices.")
  )

  # Return process
  list_to_return <- list(
    graph = generate_graph(
      ordered_distances = spanning_tree,
      number_edges = modified_output$max_value,
      vertices_names = vertices_names
    ),
    number_edges = modified_output$max_value,
    trace = modified_output$trace,
    ordered_distances = spanning_tree
  )
  class(list_to_return) <- "stad"

  return(list_to_return)
}


# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
#' Modify distance matrix
#'
#' Internal function that transforms the distance matrix into a dataframe. Additionaly, this function provides the array vertices_names for the graph generation.
#'
#' @param distance_matrix dist object containing a distance matrix.
#' @importFrom rlang .data
#'
#' @return Returns a \code{list} with the following items:
#' \itemize{
#' \item distance_matrix_df
#' \item vertices_names
#' }
modify_distance_matrix <- function(distance_matrix) {

  # Transformation of distance matrix dist into a matrix
  if ( class(distance_matrix) == "dist" ) {
    distance_matrix = as.matrix(distance_matrix)
    distance_matrix = distance_matrix / max(distance_matrix)
  } else if ( class(distance_matrix) != "matrix") {
    stop("Distance matrix must be a dist or matrix")
  }

  # Names of vertices
  vertices_names = 1:nrow(as.matrix(distance_matrix))

  # Distance matrix to Data frame structure with 'from' 'to' 'value' structure
  distance_matrix_df = reshape2::melt(as.matrix(distance_matrix), varnames = c("from", "to"))
  # Upper part of the square matrix
  distance_matrix_df = distance_matrix_df %>% dplyr::filter(.data$from < .data$to)

  # List object with the computed elements
  return(list(
    distance_matrix = distance_matrix,
    distance_matrix_df = distance_matrix_df,
    vertices_names = vertices_names
  ))
}


# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
#' Generate graph
#'
#' Internal function to generate graph from ordered_distances data.frame using the number of edges
#'
#' @param ordered_distances data.frame with the information from distance matrix and MST.
#' @param number_edges numeric value that represents the number of edges needed to generate the graph.
#' @param vertices_names numeric array with names of nodes (\code{colnames(distance_matrix)})
#' @importFrom rlang .data
#'
#' @return Returns a \code{igraph} object with the defined graph.
generate_graph <- function(ordered_distances = NULL,
                           number_edges = NULL,
                           vertices_names = NULL){

  # Filter the distances with the number of edges selected
  filtered_ordered_distances <- ordered_distances %>%
    dplyr::filter( order <= number_edges | .data$mst == 1 )

  # Define a new graph
  graph = igraph::graph_from_data_frame(filtered_ordered_distances,
                                        directed = FALSE,
                                        vertices = vertices_names)

  # Add weight
  igraph::E(graph)$weight = filtered_ordered_distances$value + 1e-8

  return(graph)
}

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
#' Modify output from optim_sa
#'
#' Transform the output from the intermediate results and save it as a \code{data.frame}.
#'
#' @param output String with the intermediate results from optim_sa.
#'
#' @return Returns a \code{list} object with following items:
#' \itemize{
#' \item max_value
#' \item trace
#' }
modify_output <- function(output) {

  # Internal function that removes the string [1]
  remove_brackets_num <- function  (string) { gsub("\\[1\\]", "", string) }

  # Detect the index where the optim_sa output (final result) starts
  par_idx <- which(output == "$par")
  par_value <- as.numeric( remove_brackets_num( output[which(output == "$par") +1 ] ) )

  # Adapt the trace from intermediate results
  intermediate <- unlist(strsplit( remove_brackets_num( output[1:par_idx-1] ), " "))
  intermediate_matrix <- as.data.frame(matrix( as.numeric(intermediate[which(intermediate != "")]) , ncol = 2, byrow = TRUE))
  colnames(intermediate_matrix) <- c("iteration", "correlation")

  return(list(max_value = par_value, trace = intermediate_matrix))
}



