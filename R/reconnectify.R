
Reconnectify <- function(g) {
  full.g <- graph_from_adjacency_matrix(xy.dist,
                                        weighted=TRUE)
  mst <- minimum.spanning.tree(full.g,
                               weights = E(full.g)$weight)

  sorted_edges <- order(E(mst)$weight)
  edge_list <- get.edgelist(mst) |> as.data.frame()

  comp <- components(knng)
  ncomp <- comp$no
  for(i in sorted_edges) {
    #print(E(mst)[i]$weight)
    g.new <- add_edges(knng, edge_list[i,],
                       weight=E(mst)[i]$weight)
    print(sum(is.na(E(g.new)$weight)))
    if(components(g.new)$no < ncomp) {
      ncomp <- ncomp - 1
      knng <- g.new
    }
    if(ncomp == 1) {
      break
    }
  }

  return(knng)
}




