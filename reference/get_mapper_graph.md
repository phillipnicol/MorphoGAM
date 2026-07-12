# Mapper Algorithm

See the following article for more details: Singh, Gurjeet, Facundo
Mémoli, and Gunnar E. Carlsson. "Topological methods for the analysis of
high dimensional data sets and 3d object recognition." SPBG. 2007.

## Usage

``` r
get_mapper_graph(
  graph,
  filter_values,
  num_intervals = 10,
  percent_overlap = 33
)
```

## Arguments

- graph:

  An igraph object (n vertices).

- filter_values:

  A numeric vector (n) of filter values, specifying how to order the
  vertices of the graph and thus group them into level sets.

- num_intervals:

  An integer specifying the resolution of the Mapper algorithm (the
  number of level sets).

- percent_overlap:

  A numeric between 0 and 100 specifying the gain of Mapper algorithm
  (the percentage that level sets should overlap).

## Value

igraph object containing the Mapper graph

## Examples

``` r
my_data <- makeSwissRoll()
#> Error in makeSwissRoll(): could not find function "makeSwissRoll"
my_graph <- dimRed:::makeKNNgraph(x = my_data, k = 5, eps = 0)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'nrow': object 'my_data' not found
mapper_graph <- get_mapper_graph(my_graph, my_data[, 1])
#> Error in get_mapper_graph(my_graph, my_data[, 1]): could not find function "get_mapper_graph"
get_mapper_plot(my_data, mapper_graph)
#> Error in get_mapper_plot(my_data, mapper_graph): could not find function "get_mapper_plot"
```
