# Predict Matches Based on a Given Record Linkage Model

Predicts matches between records in two datasets based on a given record
linkage model, using the maximum entropy classification (MEC) algorithm
(see [Lee et al.
(2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm)).

## Usage

``` r
# S3 method for class 'rec_lin_model'
predict(
  object,
  newdata_A,
  newdata_B,
  duplicates_in_A = FALSE,
  set_construction = c("size", "flr", "mmr"),
  fixed_method = "Newton",
  target_rate = 0.03,
  tol = 0.005,
  max_iter = 50,
  data_type = c("data.frame", "data.table", "matrix"),
  true_matches = NULL,
  ...
)
```

## Arguments

- object:

  A `rec_lin_model` object from the `train_rec_lin` or
  `custom_rec_lin_model` functions.

- newdata_A:

  A duplicate-free `data.frame` or `data.table`.

- newdata_B:

  A duplicate-free `data.frame` or `data.table`.

- duplicates_in_A:

  Logical indicating whether to allow `A` to have duplicate records.

- set_construction:

  A method for constructing the predicted set of matches (`"size"`,
  `"flr"` or `"mmr"`).

- fixed_method:

  A method for solving fixed-point equations using the
  [FixedPoint](https://rdrr.io/pkg/FixedPoint/man/FixedPoint.html)
  function.

- target_rate:

  A target false link rate (FLR) or missing match rate (MMR) (used only
  if `set_construction == "flr"` or `set_construction == "mmr"`).

- tol:

  Error tolerance in the bisection procedure (used only if
  `set_construction == "flr"` or `set_construction == "mmr"`).

- max_iter:

  A maximum number of iterations for the bisection procedure (used only
  if `set_construction == "flr"` or `set_construction == "mmr"`).

- data_type:

  Data type for predictions with a custom ML model (`"data.frame"`,
  `"data.table"` or `"matrix"`; used only if `object` is from the
  `custom_rec_lin_model` function).

- true_matches:

  A `data.frame` or `data.table` indicating true matches.

- ...:

  Additional controls passed to the `predict` function for custom ML
  model (used only if the `object` is from the `custom_rec_lin_model`
  function).

## Value

Returns a list containing:  

- `M_est` – a `data.table` with predicted matches,

- `set_construction` – a method for constructing the predicted set of
  matches,

- `n_M_est` – estimated classification set size,

- `flr_est` – estimated false link rate (FLR),

- `mmr_est` – estimated missing match rate (MMR),

- `iter` – the number of iterations in the bisection procedure,

- `eval_metrics` – standard metrics for quality assessment, if
  `true_matches` is provided,

- `confusion` – confusion matrix, if `true_matches` is provided.

## Details

The `predict` function estimates the probability/density ratio between
matches and non-matches for pairs in given datasets, based on a model
obtained using the `train_rec_lin` or `custom_rec_lin_model` functions.
Then, it estimates the number of matches and returns the predicted
matches, using the maximum entropy classification (MEC) algorithm (see
[Lee et al.
(2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm)).

The `predict` function allows the construction of the predicted set of
matches using its estimated size or the bisection procedure, described
in [Lee et al.
(2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm),
based on a target False Link Rate (FLR) or missing match rate (MMR). To
use the second option, set `set_construction = "flr"` or
`set_construction = "mmr"` and specify a target error rate using the
`target_rate` argument.

By default, the function assumes that the datasets `newdata_A` and
`newdata_B` contain no duplicate records. This assumption might be
relaxed by allowing `newdata_A` to have duplicates. To do so, set
`duplicates_in_A = TRUE`.

## References

Lee, D., Zhang, L.-C. and Kim, J. K. (2022). Maximum entropy
classification for record linkage. Survey Methodology, Statistics
Canada, Catalogue No. 12-001-X, Vol. 48, No. 1.

Vo, T. H., Chauvet, G., Happe, A., Oger, E., Paquelet, S., and Garès, V.
(2023). Extending the Fellegi-Sunter record linkage model for mixed-type
data with application to the French national health data system.
Computational Statistics & Data Analysis, 179, 107656.

Sugiyama, M., Suzuki, T., Nakajima, S. et al. Direct importance
estimation for covariate shift adaptation. Ann Inst Stat Math 60,
699–746 (2008).
[doi:10.1007/s10463-008-0197-x](https://doi.org/10.1007/s10463-008-0197-x)

## Author

Adam Struzik

## Examples

``` r
df_1 <- data.frame(
  "name" = c("James", "Emma", "William", "Olivia", "Thomas",
  "Sophie", "Harry", "Amelia", "George", "Isabella"),
  "surname" = c("Smith", "Johnson", "Brown", "Taylor", "Wilson",
  "Davis", "Clark", "Harris", "Lewis", "Walker")
)
 df_2 <- data.frame(
  "name" = c("James", "Ema", "Wimliam", "Olivia", "Charlotte",
  "Henry", "Lucy", "Edward", "Alice", "Jack"),
  "surname" = c("Smith", "Johnson", "Bron", "Tailor", "Moore",
  "Evans", "Hall", "Wright", "Green", "King")
)
comparators <- list("name" = jarowinkler_complement(),
                    "surname" = jarowinkler_complement())
matches <- data.frame("a" = 1:4, "b" = 1:4)
methods <- list("name" = "continuous_nonparametric",
                "surname" = "continuous_nonparametric")
model <- train_rec_lin(A = df_1, B = df_2, matches = matches,
                       variables = c("name", "surname"),
                       comparators = comparators,
                       methods = methods)
#> Warning: There are duplicate values in 'sigma', only the unique values are used.
#> Warning: There are duplicate values in 'sigma', only the unique values are used.

df_new_1 <- data.frame(
  "name" = c("John", "Emily", "Mark", "Anna", "David"),
  "surname" = c("Smith", "Johnson", "Taylor", "Williams", "Brown")
)
df_new_2 <- data.frame(
  "name" = c("John", "Emely", "Mark", "Michael"),
  "surname" = c("Smitth", "Johnson", "Tailor", "Henders")
)
predict(model, df_new_1, df_new_2)
#> The algorithm predicted 3 matches.
#> The first 3 predicted matches are:
#>        a     b ratio / 1000
#>    <num> <num>        <num>
#> 1:     2     2   0.04572852
#> 2:     3     3   0.02813814
#> 3:     1     1   0.02560156
#> ========================================================
#> The construction of the classification set was based on estimates of its size.
#> Estimated false link rate (FLR): 15.3038 %.
#> Estimated missing match rate (MMR): 15.3038 %.
```
