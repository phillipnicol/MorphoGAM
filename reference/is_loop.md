# Determines if a k-nearest neighbors graph of a dataset forms a loop.

Loops are detected as large cycles in the output of the Mapper
algorithm.

## Usage

``` r
is_loop(xy, knn_graph)
```

## Arguments

- xy:

  A numeric matrix or data frame with two columns representing the x and
  y coordinates of the data points.

- knn_graph:

  A k-nearest neighbors graph of \`xy\`.

## Value

Whether the k-nearest neighbors graph is a loop.

## Examples

``` r
my_data <- makeSwissRoll()
#> Error in makeSwissRoll(): could not find function "makeSwissRoll"
my_graph <- dimRed:::makeKNNgraph(x = my_data, k = 5, eps = 0)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'nrow': object 'my_data' not found
testthat::expect_false(is_loop(my_data, my_graph))
#> Error in is_loop(my_data, my_graph): could not find function "is_loop"
```
