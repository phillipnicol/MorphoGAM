get_clusters_in_level <- function(data, is_in_level, eps, min_pts) {
  if (!any(is_in_level)) { return(list()) }

  level_data <- data[is_in_level, ]
  level_dbscan <- dbscan::dbscan(level_data, eps = eps, minPts = min_pts)
  cluster_indices <- level_dbscan$cluster
  num_clusters <- max(cluster_indices)

  clusters <- lapply(seq_len(num_clusters), function(cluster_index) {
    which(is_in_level)[cluster_indices == cluster_index]
  })
  return(clusters)
}

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
#' @param data A numeric matrix (n x 2) representing spatial coordinates.
#' @param filter_values A numeric vector (n) of the filter values on the data.
#' @param num_intervals An integer specifying the resolution of the Mapper
#'   algorithm (how many level sets to cluster)
#' @param percent_overlap A numeric between 0 and 100 specifying the gain of
#'   Mapper algorithm (what percent should level sets overlap)
#' @param eps A numeric specifying the epsilon radius that defines neighboring
#'   points in DBSCAN.
#' @param min_pts An integer specifying the minimum cluster size for DBSCAN.
#'
#' @return igraph object containing the Mapper graph
#'
#' @examples
#' my_data <- makeSwissRoll()
#' scaled_data <- scale(my_data)
#' pc <- prcomp(scaled_data, center = FALSE, scale = FALSE)
#' mapper_graph <- get_mapper_graph(scaled_data, pc$x[, 1], num_intervals = 7,
#'                                  eps = 0.2, min_pts = 5)
#' get_mapper_plot(scaled_data, mapper_graph)
get_mapper_graph <- function(
    data,
    filter_values,
    num_intervals = 10,
    percent_overlap = 30,
    eps = 0.1,
    min_pts = 50
) {
  if (percent_overlap >= 50) {
    warning("Overlap should be under 50%.")
  }

  filter_min <- min(filter_values)
  filter_max <- max(filter_values)
  overlap <- percent_overlap / 100
  interval_span <- (num_intervals - (num_intervals - 1) * overlap)
  interval_length <- (filter_max - filter_min) / interval_span
  step_size <- interval_length * (1 - overlap)

  # Iterate through each level set and cluster the data.
  mapper_graph <- igraph::make_empty_graph(directed = FALSE)
  for (level in 1:num_intervals) {
    lo <- filter_min + (level - 1) * step_size
    hi <- lo + interval_length
    is_in_level <- (lo <= filter_values) & (filter_values <= hi)
    level_clusters <- get_clusters_in_level(data, is_in_level, eps, min_pts)

    # Add each cluster as a vertex in the graph.
    for (cluster in level_clusters) {
      attributes <- list(level=level, points=list(cluster))
      mapper_graph <- igraph::add_vertices(mapper_graph, 1, attr=attributes)
    }

    # Add edges to vertices from the previous level set.
    edges <- get_mapper_edges(igraph::V(mapper_graph), level - 1, level)
    mapper_graph <- igraph::add_edges(mapper_graph, edges)
  }
  return(mapper_graph)
}

get_mapper_plot <- function(data, mapper_graph) {
  points <- igraph::V(mapper_graph)$points
  x_coords <- lapply(points, function(x) mean(data[x, 1]))
  y_coords <- lapply(points, function(x) mean(data[x, 2]))
  coordinates <- cbind(unlist(x_coords), unlist(y_coords))
  plot(data[, 1], data[, 2])
  plot(mapper_graph, layout=coordinates, add = TRUE, rescale = FALSE,
       vertex.label = NA, vertex.size = 15, vertex.color = "orange",
       edge.color = "red")
}
