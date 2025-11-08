# Absolute Distance Comparison Function

Creates a function that calculates the absolute distance between two
values.

## Usage

``` r
abs_distance()
```

## Value

Returns a function taking two arguments, `x` and `y`, and returning
their absolute difference.

## Author

Adam Struzik

## Examples

``` r
cmp <- abs_distance()
cmp(1, 5) # returns 4
#> [1] 4
```
