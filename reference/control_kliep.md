# Controls for the `kliep` Function

Controls for the
[kliep](https://thomvolker.github.io/densityratio/reference/kliep.html)
function used in the package.

## Usage

``` r
control_kliep(scale = NULL, progressbar = FALSE, nfold = 2, ...)
```

## Arguments

- scale:

  `"numerator"`, `"denominator"` or `NULL`, indicating whether to
  standardize each numeric variable according to the numerator means and
  standard deviations, the denominator means and standard deviations, or
  apply no standardization at all.

- progressbar:

  Logical indicating whether or not to display a progressbar.

- nfold:

  Number of cross-validation folds used in order to calculate the
  optimal sigma value (default is 2-fold cv).

- ...:

  Additional arguments.

## Value

Returns a list with parameters.

## Author

Adam Struzik
