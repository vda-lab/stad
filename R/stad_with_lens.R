
# ----------------------
# Table of contents
# ----------------------
# stad_with_lens
# bivariate_split
# univariate_intra
# ----------------------

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
#' STAD with lens
#'
#' The function computes STAD with lens using a given distance matrix and filter values.
#'
#' @param distance_matrix dist object containing a distance matrix.
#' @param filter_values numeric array or two-dimensional data.frame with the lenses to use in STAD.
#' @param num_intervals number or array defining the number of intervals to consider in the filters.
#' @param metric string or array defining the metrics supported ("polar" or "euclidean").
#' @param penalty numeric penalty value used during the optimization loop. It limits the number of iterations. Default = 0.
#' @param random_factor integer factor that controls the variability of the next iteration on the optimization algorithm. Default = 10000. Higher values generates more diverse exploration.
#' @param iterations_inner_loop integer number of evaluations for each iteration. The number of iterations allows knowing the close values around but increase the total number of iterations.
#' @param two_mst boolean indicating if the MST is build using the two-step MST. It only applies when \code{filter_values} are defined.
#' @importFrom rlang .data
#'
#' @return Returns a \code{stad} class (or \code{list}) with the following items:
#' \itemize{
#' \item graph
#' \item number_edges
#' \item trace
#' \item ordered_distances
#' }
stad_with_lens <- function( distance_matrix,
                            filter_values = NULL,
                            num_intervals = NULL,
                            metric = NULL,
                            penalty = 0,
                            random_factor = 10000,
                            iterations_inner_loop = 10,
                            two_mst = FALSE){


  # Distance matrix transformation
  modified_distance_matrix <- modify_distance_matrix(distance_matrix)
  vertices_names <- modified_distance_matrix$vertices_names
  distance_matrix <- modified_distance_matrix$distance_matrix
  distance_matrix_df <- modified_distance_matrix$distance_matrix_df

  # ---------- Filter distance matrix ----------
  # Each edge is the union of two points (from, to). From and To represents two points in the high-dimensional space.
  # We use the the filter function to limit the edge evaluation, ex. if the edge from-to is to far in filter function space,
  # the distance is not evaluated in STAD
  if (class(filter_values) == "data.frame") {

    # From coordinates into filter_values space (categorized)
    from  <- data.frame(list(
      from = as.numeric(colnames(distance_matrix)),
      from_x = as.numeric(as.factor(cut(filter_values$x, num_intervals[1], labels = FALSE, include.lowest = TRUE))),
      from_y = as.numeric(as.factor(cut(filter_values$y, num_intervals[2], labels = FALSE, include.lowest = TRUE)))
    ))
    # To coordinates into filter_values space (categorized)
    to <- from %>% dplyr::rename(to = from, to_x = .data$from_x, to_y = .data$from_y)

    # Define default metric
    if( is.null(metric) ){
      metric <- c("euclidean", "euclidean")
    } else if ( length(metric) == 1 ) {
      stop("'metric' parameter has dimension 1 but 'filter_values' uses 2 dimensions.")
    }

    # Merge
    # We filter all values outside bivariate_split result
    distance_matrix_df  <- distance_matrix_df %>%
      dplyr::inner_join(from, by  = "from") %>%
      dplyr::inner_join(to, by = "to") %>%
      dplyr::left_join(bivariate_split(from, metric = metric), by = c("from_x", "from_y", "to_x", "to_y")) %>%
      dplyr::mutate(intra = ifelse(.data$from_x == .data$to_x & .data$from_y == .data$to_y, paste0(.data$from_x , .data$from_y), NA),
             inter = ifelse( !is.na(.data$inter), paste0(.data$from_x , .data$from_y), NA))
  } else if (class(filter_values) == "numeric") {
    # From
    from <- data.frame(list(
      from = as.numeric(colnames(distance_matrix)),
      from_x = as.numeric(as.factor(cut(filter_values, num_intervals, labels = FALSE, include.lowest = TRUE)))  )
    )
    # To
    to <- from %>% dplyr::rename(to = from, to_x = .data$from_x)

    # Define default metric
    if( is.null(metric) ){
      metric <- "euclidean"
    } else if ( length(metric) == 2 ) {
      stop("'metric' parameter has dimension 2 but 'filter_values' just 1 dimension.")
    }

    # Range of values for univariate intra function
    range_groups <- range(c(from$from_x, to$to_x))

    # Merge
    distance_matrix_df <- distance_matrix_df %>%
      dplyr::inner_join(from, by  = "from") %>%
      dplyr::inner_join(to, by = "to") %>%
      dplyr::mutate(
        inter = ifelse(.data$from_x == .data$to_x, .data$from_x, NA),
        intra = univariate_intra(
          .data$from_x,
          .data$to_x,
          range_groups,
          metric
        )
      )
  } else {
    stop("'filter_values' format is not supported by STAD. It only supports 'data.frame' and 'numeric' classes.")
  }

  # Define IDs to be assigned to NA's. Values out of lens scope
  ids_NA <- distance_matrix_df %>%
    dplyr::filter(is.na(.data$inter) & is.na(.data$intra)) %>%
    dplyr::select(from, to)

  # Assign NA's to values out of lens scope
  # New distance matrix (matrix fromat)
  distance_matrix[as.matrix(ids_NA)] <- NA

  # New distance matrix (data.frame format)
  # Weight inter is increased to 2 + value to prioritize intra
  distance_matrix_df <-  distance_matrix_df %>%
    dplyr::filter(!is.na(.data$inter) | !is.na(.data$intra)) %>%
    dplyr::mutate(
      # We apply this correction to the inter-edges.
      # This is the new weight that we will use in the Step 1.
      weight_intra = ifelse(!is.na(.data$inter), .data$value, 2 + .data$value)
    )
  # ---------- End filter distance matrix ----------

  if(two_mst == TRUE){
    # ---------- Graph transfromation and MST computation ----------
    # Step 1: Intra-edges
    # Distance matrix is transformed into a graph (igraph)
    distance_matrix_graph_intra <- igraph::graph_from_data_frame(distance_matrix_df, directed = FALSE, vertices = vertices_names)

    # Computing MST using within distance and 2 + between distance to prioritize intra-edges
    # Weight for the graph using inter values (inter = 2 + values and intra = values)
    igraph::E(distance_matrix_graph_intra)$weight <- distance_matrix_df$weight_intra
    # Computing the MST with the defined weight
    mst_graph_intra <- igraph::mst(distance_matrix_graph_intra, weights = distance_matrix_df$weight_intra)

    # Cluster random walk
    walktrap_intra <- igraph::cluster_walktrap(
      weights = 1 / (1 + distance_matrix_df$weight_intra),
      mst_graph_intra,
      steps = 20
    )
    # Define cluster_intra_from as a data.frame with the list of vertices and the clusters
    cluster_intra_from <- as.data.frame(cbind(
      vertices_names,
      as.matrix(igraph::membership(walktrap_intra))
    ))
    names(cluster_intra_from) <- c("from", "group_from")
    cluster_intra_to <- cluster_intra_from; names(cluster_intra_to) = c("to", "group_to")

    # Step 2: Inter-edges
    # Distance Matrix With the Edges from MST-intra
    # Adding clusters to mst_intra
    distance_matrix_mst_intra <- igraph::as_data_frame(mst_graph_intra) %>%
      dplyr::mutate(from = as.numeric(from),
                    to = as.numeric(to) ) %>%
      dplyr::inner_join(cluster_intra_from, by = "from") %>%
      dplyr::inner_join(cluster_intra_to, by = "to") %>%
      dplyr::full_join(
        # We add all intra connections into the MST. In the previous steps,
        # we penalize the inter distances. In  the second step, we will
        # evaluate inter connections (without penalty).
        # We extract the inter-edges connections from the original distance matrix.
        distance_matrix_df %>% dplyr::filter(!is.na(.data$intra)),
        by = c("from", "to", "value", "from_x", "to_x", "inter", "intra", "weight_intra")
      ) %>%
      dplyr::mutate(weight_inter = ifelse(
        # If the nodes in the edge belong to the same cluster AND this edge is intra-edge
        # then we will preserve this value.
        # Remaining values will be evaluated right after. These values are in the range
        # weight = value + 2 to keep selected edges (distances = values).
        .data$group_from == .data$group_to & !is.na(.data$inter),
        .data$value,
        .data$value + 2
      ))

    # MST to build the definitive inter-edges but keeping the selected intra-edges
    distance_matrix_graph_inter <- igraph::graph_from_data_frame(
      distance_matrix_mst_intra,
      directed = FALSE,
      vertices = vertices_names
    )
    mst_graph_inter <- igraph::mst(
      distance_matrix_graph_inter,
      weights = distance_matrix_mst_intra$weight_inter
    )
  } else if (two_mst == FALSE) {

    # Distance matrix is transformed into a graph (igraph)
    distance_matrix_graph_intra <- igraph::graph_from_data_frame(distance_matrix_df, directed = FALSE, vertices = vertices_names)
    # Computing the MST with the defined weight
    mst_graph_inter <- igraph::mst(
      distance_matrix_graph_intra,
      weights = distance_matrix_df$value
    )

  }

  # ---------- End graph transfromation and MST computation ----------

  # STAD iteration
  list_to_return <- stad_evaluation(distance_matrix = distance_matrix,
                                    distance_matrix_df = distance_matrix_df,
                                    vertices_names = vertices_names,
                                    mst_graph = mst_graph_inter,
                                    penalty = penalty,
                                    random_factor = random_factor,
                                    iterations_inner_loop = iterations_inner_loop)

  return(list_to_return)
}


# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
#' Bivariate split
#'
#' Filter values connections for two-dimensional lens
#'
#' @param split data.frame with the connections (categories) of two-dimensional lenses
#' @param x string variable name for the first dimension of the lens. Default 'from_x'.
#' @param y string variable name for the second dimension of the lens. Default 'from_y'.
#' @param metric string or array defining the metrics supported ("polar" or "euclidean").
#' @importFrom rlang .data
#'
#' @return Returns a \code{data.frame} with the preserved connections.
bivariate_split <-  function (split, x = "from_x", y = "from_y", metric = NULL) {

  # Create grouping. We define a unique set of combinations between the filter values x and y
  connections <- split %>% dplyr::select_(x, y) %>% dplyr::rename_( "x" = x, "y" = y)
  connections_group <- as.matrix( unique ( connections ) )

  # Compute distance matrix. The output generates a distance matrix with the selected metric
  distance_matrix <- as.matrix(custom_distance(connections_group, metric))

  # Names of vertices
  vertices_names <- 1:nrow(as.matrix(distance_matrix))
  # rename distance_matrix
  colnames(distance_matrix) <- vertices_names
  rownames(distance_matrix) <- vertices_names

  # DF of distance matrix
  distance_matrix_df <- reshape2::melt(distance_matrix, varnames = c("from", "to"))

  # Distance matrix is transformed into a graph (igraph)
  distance_matrix_graph <- igraph::graph_from_data_frame(distance_matrix_df, directed = FALSE, vertices = vertices_names)

  # Computing Minimum Spanning Tree from igraph package (MST) to detect which combinations of filter_values are closest.
  # The MST allows a connection between all elements
  mst_matrix <- igraph::as_data_frame(
    igraph::mst(distance_matrix_graph, weights = distance_matrix_df$value)
  ) %>%
    dplyr::mutate(from = as.numeric(from),
                  to = as.numeric(to),
                  mst = 1)

  # Complentary connection --> ex. connection 1,2 and 2,1
  mst_matrix_inverse <- mst_matrix %>%
    dplyr::mutate(from2 = from,
                  from = to,
                  to = .data$from2) %>%
    dplyr::select(-.data$from2)

  # Merge. All diagonal connections are included in the final graph.
  spanning_tree <- dplyr::full_join(distance_matrix_df,
                                    dplyr::bind_rows(mst_matrix, mst_matrix_inverse),
                                    by = c("from", "to", "value") ) %>%
    dplyr::mutate(
      mst = ifelse(is.na(.data$mst), 0, .data$mst),
      value_mst = ifelse(.data$mst == 1, sqrt(2*.data$value), NA)) %>%
    dplyr::filter(from != to) %>%
    dplyr::group_by(from) %>%
    dplyr::mutate(value_mst = min(.data$value_mst, na.rm = TRUE) ) %>%
    dplyr::ungroup() %>%
    dplyr::filter( (.data$value <= .data$value_mst | .data$mst == 1) ) # & from < to

  # Adding previous information
  # from
  from <- as.data.frame(connections_group) %>%
    dplyr::mutate(from = dplyr::row_number()) %>%
    dplyr::rename(from_x = x, from_y = y)
  # to
  to <-  from %>% dplyr::rename(to = from, to_x = .data$from_x, to_y = .data$from_y)
  # merge
  spanning_tree  = spanning_tree %>%
    dplyr::inner_join(from, by  = "from") %>%
    dplyr::inner_join(to, by = "to") %>%
    dplyr::mutate(inter = 1) %>%
    dplyr::select(.data$from_x, .data$from_y, .data$to_x, .data$to_y, .data$inter)

  return(spanning_tree)
}

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
#' Custom distance
#'
#' Internal distance matrix for bivariate_split. Uses polar or euclidean metric. Returns a distance matrix as sum of the two independents dimensions.
#'
#' @param x array variable. Dimension of the lens.
#' @param metric string or array defining the metrics supported ("polar" or "euclidean").
#'
#' @return Returns a \code{dist} object with the distance of the filter.
custom_distance <- function (x, metric) {

  # Scale values between 0 and 1
  range01 <- function(x){
    range <- range(x)
    return((x - range[1]) / (range[2] - range[1]))
  }

  # Custom polar distance computation
  distance_polar <- function(v1, v2) {
    # Degree to radian
    deg2rad <- function(x) { (2 * x * pi) }
    return(
      sqrt(2  - 2 * cos( deg2rad( v1 - v2 ) )) / 2
    )
  }

  # Standardization
  x <- apply(x, 2 , range01)
  x_1 <- x[,1]
  x_2 <- x[,2]

  if(metric[1] == "euclidean" & metric[2] == "euclidean") {
    x_1_dist = as.matrix(stats::dist(as.matrix(x_1), method = "euclidean"))
    x_2_dist = as.matrix(stats::dist(as.matrix(x_2), method = "euclidean"))
  } else if(metric[1] == "polar" & metric[2] == "euclidean") {
    x_1_dist = as.matrix(usedist::dist_make(as.matrix(x_1), distance_polar))
    x_2_dist = as.matrix(stats::dist(as.matrix(x_2), method = "euclidean"))
  } else if(metric[1] == "euclidean" & metric[2] == "polar") {
    x_1_dist = as.matrix(stats::dist(as.matrix(x_1), method = "euclidean"))
    x_2_dist = as.matrix(usedist::dist_make(as.matrix(x_2), distance_polar))
  } else if(metric[1] == "polar" & metric[2] == "polar") {
    x_1_dist = as.matrix(usedist::dist_make(as.matrix(x_1), distance_polar))
    x_2_dist = as.matrix(usedist::dist_make(as.matrix(x_2), distance_polar))
  }

  return(stats::as.dist( sqrt(x_1_dist**2 + x_2_dist**2) ))
}

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
#' Univariate intra-edge
#'
#' Definition of intra-edge for one-dimensional lens
#'
#' @param from_x numeric; categroy of the group (interval) for from element in the edge.
#' @param to_x numeric; categroy of the group (interval) for to element in the edge.
#' @param range array. Numeric range of the lens
#' @param metric string or array defining the metrics supported ("polar" or "euclidean").
#'
#' @return Returns a logical element.
univariate_intra <- function(from_x, to_x, range, metric){

  if (metric == "euclidean" ) {
    intra <- ifelse(from_x + 1 == to_x | from_x == to_x + 1, from_x, NA)
  } else if(metric == "polar") {
    intra <- ifelse(
      from_x + 1 == to_x |
        from_x == to_x + 1 |
        (from_x == range[1] & to_x == range[2]) |
        (to_x == range[1] & from_x == range[2]), from_x, NA)
  }

  return(intra)

}
