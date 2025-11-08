# Create a Custom Record Linkage Model

Creates a supervised record linkage model using a custom machine
learning (ML) classifier.

## Usage

``` r
custom_rec_lin_model(ml_model, vectors)
```

## Arguments

- ml_model:

  A trained ML model that predicts the probability of a match based on
  comparison vectors.

- vectors:

  An object of class `comparison_vectors` (a result of the
  `comparison_vectors` function), used for training the `ml_model`.

## Value

Returns a list containing:  

- `b_vars` – here `NULL`,

- `cpar_vars` – here `NULL`,

- `cnonpar_vars` – here `NULL`,

- `b_params` – here `NULL`,

- `cpar_params` – here `NULL`,

- `cnonpar_params` – here `NULL`,

- `ratio_kliep` – here `NULL`,

- `ratio_kliep_list` – here `NULL`,

- `ml_model` – ML model used for creating the record linkage model,

- `pi_est` – a prior probability of matching,

- `match_prop` – proportion of matches in the smaller dataset,

- `variables` – a character vector of key variables used for comparison,

- `comparators` – a list of functions used to compare pairs of records,

- `methods` – here `NULL`,

- `prob_ratio` – here `"2"`.

## Details

The `custom_rec_lin_model` function creates a custom record linkage
model, based on known matches and non-matches (which might later serve
as a classifier for pairs outside training data). The procedure of
creating a custom model based on training data is as follows.

1.  Use the `comparison_vectors` function to compare pairs of records.

2.  Train a machine learning classifier using the `Omega` element of the
    output of the `comparison_vectors` function. The classifier should
    predict the probability of matching based on a given vector.

3.  Use the `custom_rec_lin_model` function with appropriate arguments.

## Author

Adam Struzik

## Examples

``` r
if (requireNamespace("xgboost", quietly = TRUE)) {
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
  vectors <- comparison_vectors(A = df_1, B = df_2, variables = c("name", "surname"),
                               comparators = comparators, matches = matches)
  train_data <- xgboost::xgb.DMatrix(
    data = as.matrix(vectors$Omega[, c("gamma_name", "gamma_surname")]),
    label = vectors$Omega$match
  )
  params <- list(objective = "binary:logistic",
                 eval_metric = "logloss")
  model_xgb <- xgboost::xgboost(data = train_data, params = params,
                                nrounds = 100, verbose = 0)
  custom_xgb_model <- custom_rec_lin_model(model_xgb, vectors)
  custom_xgb_model
}
#> Record linkage model based on the following variables: name, surname.
#> A custom ML model was used.
#> The prior probability of matching is 0.04.
#> ========================================================
#> Probability/density ratio type: 2.
```
