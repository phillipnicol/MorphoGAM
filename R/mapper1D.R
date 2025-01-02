# Returns the components in a graph's level set.
get_clusters_in_level <- function(graph, is_in_level) {
  if (!any(is_in_level)) { return(list()) }

  level_vertices <- igraph::V(graph)[is_in_level]
  level_subgraph <- igraph::induced_subgraph(graph, level_vertices)
  level_components <- igraph::components(level_subgraph)

  clusters <- lapply(seq_along(level_components$csize), function(component) {
    is_in_component <- (level_components$membership == component)
    names(igraph::V(level_subgraph)[is_in_component])
  })
  return(clusters)
}

# Returns the pairs of vertices from different level sets that share a point.
get_mapper_edges <- function(vertices, level1, level2) {
  clusters_overlap <- function(v1, v2) {
    cluster1 <- unlist(vertices[v1]$points)
    cluster2 <- unlist(vertices[v2]$points)
    length(intersect(cluster1, cluster2)) > 0
  }
  edges <- lapply(vertices[vertices$level == level1], function(v1) {
    lapply(vertices[vertices$level == level2], function(v2) {
      if (clusters_overlap(v1, v2)) { c(v1, v2) }
    })
  })
  unlist(edges)
}

#' Mapper Algorithm
#'
#' See the following article for more details:
#' Singh, Gurjeet, Facundo Mémoli, and Gunnar E. Carlsson. "Topological methods
#' for the analysis of high dimensional data sets and 3d object recognition."
#' SPBG. 2007.
#'
#' @param graph An igraph object (n vertices).
#' @param filter_values A numeric vector (n) of filter values, specifying how to
#'   order the vertices of the graph and thus group them into level sets.
#' @param num_intervals An integer specifying the resolution of the Mapper
#'   algorithm (the number of level sets).
#' @param percent_overlap A numeric between 0 and 100 specifying the gain of
#'   Mapper algorithm (the percentage that level sets should overlap).
#'
#' @return igraph object containing the Mapper graph
#'
#' @examples
#' my_data <- makeSwissRoll()
#' my_graph <- dimRed:::makeKNNgraph(x = my_data, k = 5, eps = 0)
#' mapper_graph <- get_mapper_graph(my_graph, my_data[, 1])
#' get_mapper_plot(my_data, mapper_graph)
get_mapper_graph <- function(
    graph,
    filter_values,
    num_intervals = 10,
    percent_overlap = 33
) {
  if (percent_overlap >= 50) {
    warning("Overlap should be under 50%.")
  }
  # Vertex names are necessary to identify subgraph components.
  vertex_names <- names(igraph::V(graph))
  if (is.null(vertex_names)) {
    igraph::V(graph)$name <- as.character(igraph::V(graph))
  }

  filter_min <- min(filter_values)
  filter_max <- max(filter_values)
  overlap <- percent_overlap / 100
  interval_scaling_factor <- (num_intervals - (num_intervals - 1) * overlap)
  interval_length <- (filter_max - filter_min) / interval_scaling_factor
  step_size <- interval_length * (1 - overlap)

  # Iterate through each level set and identify the components of its subgraph.
  mapper_graph <- igraph::make_empty_graph(directed = FALSE)
  for (level in 1:num_intervals) {
    lo <- filter_min + (level - 1) * step_size
    hi <- lo + interval_length
    is_in_level <- (lo <= filter_values) & (filter_values <= hi)
    level_clusters <- get_clusters_in_level(graph, is_in_level)

    # Add each component as a vertex in the graph.
    for (cluster in level_clusters) {
      # TODO: Consider ignoring small clusters.
      # if (length(cluster) < 5) next
      attributes <- list(level=level, points=list(cluster))
      mapper_graph <- igraph::add_vertices(mapper_graph, 1, attr=attributes)
    }

    # Add edges to vertices from the previous level set.
    edges <- get_mapper_edges(igraph::V(mapper_graph), level - 1, level)
    mapper_graph <- igraph::add_edges(mapper_graph, edges)
  }
  return(mapper_graph)
}

# Returns the length of the largest cycle in the graph.
#
# This is implemented inefficiently, by finding all simple paths between
# neighbors in the graph, but it should be fast enough for small graphs.
get_largest_cycle_length <- function(graph) {
  paths_to_neighbors <- lapply(igraph::V(graph), function(v) {
    curr_level <- igraph::V(graph)[v]$level
    # Performance optimization: Only consider neighbors from the previous level.
    prev_neighbors <- igraph::neighbors(graph, v)[level == curr_level - 1]
    if (!any(prev_neighbors)) { return(NULL) }
    igraph::all_simple_paths(graph, v, prev_neighbors)
  })
  paths_to_neighbors <- unlist(paths_to_neighbors, recursive = FALSE)
  path_lengths <- vapply(paths_to_neighbors, length, FUN.VALUE=numeric(1))
  # Ignore trivial paths between neighbors.
  if (all(path_lengths <= 2)) { return(0) }
  max(path_lengths[path_lengths > 2])
}

#' Determines if a k-nearest neighbors graph of a dataset forms a loop.
#'
#' Loops are detected as large cycles in the output of the Mapper algorithm.
#'
#' @param xy A numeric matrix or data frame with two columns representing the x
#'   and y coordinates of the data points.
#' @param knn_graph A k-nearest neighbors graph of `xy`.
#'
#' @return Whether the k-nearest neighbors graph is a loop.
#'
#' @examples
#' my_data <- makeSwissRoll()
#' my_graph <- dimRed:::makeKNNgraph(x = my_data, k = 5, eps = 0)
#' testthat::expect_false(is_loop(my_data, my_graph))
is_loop <- function(xy, knn_graph) {
  principal_components <- stats::prcomp(xy, center = TRUE, scale = TRUE)
  num_intervals <- 10
  mapper_graph <- get_mapper_graph(knn_graph, principal_components$x[, 1],
                                   num_intervals = num_intervals,
                                   percent_overlap = 33)
  cycle_length <- get_largest_cycle_length(mapper_graph)
  # True when a cycle takes up ~1/2 the range of the first principal component.
  return(cycle_length >= num_intervals)
}

# Temporary helper function to visualize the output of the Mapper algorithm.
get_mapper_plot <- function(xy, mapper_graph) {
  points <- lapply(igraph::V(mapper_graph)$points, as.numeric)
  x_coords <- lapply(points, function(p) mean(xy[p, 1]))
  y_coords <- lapply(points, function(p) mean(xy[p, 2]))
  coordinates <- cbind(unlist(x_coords), unlist(y_coords))
  plot(xy)
  plot(mapper_graph, layout=coordinates, add = TRUE, rescale = FALSE,
       vertex.label = NA, vertex.size = 10, vertex.color = "orange",
       edge.color = "red")
}
