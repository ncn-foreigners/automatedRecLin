The `automatedRecLin` Package
================

## Description

This R package is designed to perform record linkage (also known as
entity resolution) in unsupervised or supervised settings. It compares
pairs of records from two datasets using selected comparison functions
to estimate the probability or density ratio between matched and
non-matched records. Based on these estimates, it predicts a set of
matches that maximizes entropy.

## Installation

To install the development version from GitHub you can use the `pak`
package.

``` r
# install.packages("pak") # uncomment if needed
pak::pkg_install("ncn-foreigners/automatedRecLin")
```

## Basic usage

Load the package for the examples.

``` r
library(automatedRecLin)
```

### Unsupervised maximum entropy classifier for record linkage

Generate two simple datasets that contain some common records, with
typos in some cases.

``` r
df_1 <- data.frame(
  name = c("Emma", "Liam", "Olivia", "Noah", "Ava",
           "Ethan", "Sophia", "Mason", "Isabella", "James"),
  surname = c("Smith", "Johnson", "Williams", "Brown", "Jones",
              "Garcia", "Miller", "Davis", "Rodriguez", "Wilson"),
  city = c("New York", "Los Angeles", "Chicago", "Houston", "Phoenix",
           "Philadelphia", "San Antonio", "San Diego", "Dallas", "San Jose")
)
df_2 <- data.frame(
  name = c(
    "Emma", "Liam", "Olivia", "Noah",
    "Ava", "Ehtan", "Sopia", "Mson",
    "Charlotte", "Benjamin", "Amelia", "Lucas"
  ),
  surname = c(
    "Smith", "Johnson", "Williams", "Brown",
    "Jnes", "Garca", "Miler", "Dvis",
    "Martinez", "Lee", "Hernandez", "Clark"
  ),
  city = c(
    "New York", "Los Angeles", "Chicago", "Houston",
    "Phonix", "Philadelpia", "San Antnio", "San Dieg",
    "Seattle", "Miami", "Boston", "Denver"
  )
)
df_1
#>        name   surname         city
#> 1      Emma     Smith     New York
#> 2      Liam   Johnson  Los Angeles
#> 3    Olivia  Williams      Chicago
#> 4      Noah     Brown      Houston
#> 5       Ava     Jones      Phoenix
#> 6     Ethan    Garcia Philadelphia
#> 7    Sophia    Miller  San Antonio
#> 8     Mason     Davis    San Diego
#> 9  Isabella Rodriguez       Dallas
#> 10    James    Wilson     San Jose
df_2
#>         name   surname        city
#> 1       Emma     Smith    New York
#> 2       Liam   Johnson Los Angeles
#> 3     Olivia  Williams     Chicago
#> 4       Noah     Brown     Houston
#> 5        Ava      Jnes      Phonix
#> 6      Ehtan     Garca Philadelpia
#> 7      Sopia     Miler  San Antnio
#> 8       Mson      Dvis    San Dieg
#> 9  Charlotte  Martinez     Seattle
#> 10  Benjamin       Lee       Miami
#> 11    Amelia Hernandez      Boston
#> 12     Lucas     Clark      Denver
```

Specify the key variables used for record linkage. Select a comparison
function (i.e. a function to compare pairs of records) for each
variable. For example, use the `jarowinkler_complement` function from
the `automatedRecLin` package (1 - Jaro-Winkler distance). Choose a
method for estimating the probability or density ratio for each
variable. The available methods are: `"binary"`,
`"continuous_parametric"` and `"continuous_nonparametric"`.

``` r
variables <- c("name", "surname", "city")
comparators <- list(
  "name" = jarowinkler_complement(),
  "surname" = jarowinkler_complement(),
  "city" = jarowinkler_complement()
)
methods <- list(
  "name" = "continuous_parametric",
  "surname" = "continuous_parametric",
  "city" = "continuous_parametric"
)
```

Perform record linkage using the `mec` function. The output contains the
following information:

- the names of key variables,
- the number of predicted matches,
- the first 6 predicted matches (with their estimated probability or
  density ratio),
- the method for constructing the predicted set of matches (default:
  `"size"`),
- estimated false link rate (FLR),
- estimated missing match rate (MMR),
- estimated parameters for the variables using the `"binary"` or
  `"continuous_parametric"` methods.

``` r
set.seed(1)
unsup_result <- mec(A = df_1, B = df_2,
                    variables = variables,
                    comparators = comparators,
                    methods = methods)
unsup_result
#> Record linkage based on the following variables: name, surname, city.
#> ========================================================
#> The algorithm predicted 8 matches.
#> The first 6 predicted matches are:
#>        a     b ratio / 1000
#>    <num> <num>        <num>
#> 1:     6     6 1.433031e+08
#> 2:     8     8 3.198692e+07
#> 3:     7     7 9.673745e+05
#> 4:     5     5 2.813745e+04
#> 5:     1     1 3.375000e+00
#> 6:     2     2 3.375000e+00
#> ========================================================
#> The construction of the classification set was based on estimates of its size.
#> Estimated false link rate (FLR): 0.2066 %.
#> Estimated missing match rate (MMR): 0.0000 %.
#> ========================================================
#> Variables selected for the continuous parametric method: name, surname, city.
#> Estimated parameters for the continuous parametric method:
#>         variable p_0_M    alpha_M   beta_M      p_0_U  alpha_U    beta_U
#>           <char> <num>      <num>    <num>      <num>    <num>     <num>
#> 1:    gamma_name 0.625 138.462279 2199.107 0.04166667 6.516736 11.173089
#> 2: gamma_surname 0.500 120.665706 1974.530 0.03333333 4.622775  7.167261
#> 3:    gamma_city 0.500   6.512723  135.163 0.03333333 5.233194  9.313035
```

### Supervised maximimum entropy classifier for record linkage

Generate two simple training datasets that contain some common records,
with typos in some cases.

``` r
df_1_train <- data.frame(
        "name" = c("John", "Emily", "Mark", "Anna", "David"),
        "surname" = c("Smith", "Johnson", "Taylor", "Williams", "Brown")
)
df_2_train <- data.frame(
        "name" = c("John", "Emely", "Marc", "Michael"),
        "surname" = c("Smith", "Jonson", "Tailor", "Henderson")
)
df_1_train
#>    name  surname
#> 1  John    Smith
#> 2 Emily  Johnson
#> 3  Mark   Taylor
#> 4  Anna Williams
#> 5 David    Brown
df_2_train
#>      name   surname
#> 1    John     Smith
#> 2   Emely    Jonson
#> 3    Marc    Tailor
#> 4 Michael Henderson
```

Specify the key variables, select comparison functions and choose
methods for estimating the probability or density ratio. Additionally,
provide a `data.frame` indicating known matches.

``` r
variables_train <- c("name", "surname")
comparators_train <- list("name" = jarowinkler_complement(),
                          "surname" = jarowinkler_complement())
methods_train <- list("name" = "continuous_nonparametric",
                      "surname" = "continuous_nonparametric")
matches_train <- data.frame("a" = 1:3, "b" = 1:3)
```

Train a record linkage model using the `train_rec_lin` function.

``` r
model <- train_rec_lin(A = df_1_train, B = df_2_train,
                       matches = matches_train,
                       variables = variables_train,
                       comparators = comparators_train,
                       methods = methods_train)
model
#> Record linkage model based on the following variables: name, surname.
#> The prior probability of matching is 0.15.
#> ========================================================
#> Variables selected for the continuous nonparametric method: name, surname.
```

Generate two new datasets for record linkage prediction.

``` r
df_1_new <- data.frame(
  "name" = c("Jame", "Lia", "Tomas", "Matthew", "Andrew"),
  "surname" = c("Wilsen", "Thomsson", "Davis", "Robinson", "Scott")
)
df_2_new <- data.frame(
  "name" = c("James", "Leah", "Thomas", "Mathew", "Andrew", "Sophie"),
  "surname" = c("Wilson", "Thompson", "Davies", "Robins", "Scots", "Clarks")
)
df_1_new
#>      name  surname
#> 1    Jame   Wilsen
#> 2     Lia Thomsson
#> 3   Tomas    Davis
#> 4 Matthew Robinson
#> 5  Andrew    Scott
df_2_new
#>     name  surname
#> 1  James   Wilson
#> 2   Leah Thompson
#> 3 Thomas   Davies
#> 4 Mathew   Robins
#> 5 Andrew    Scots
#> 6 Sophie   Clarks
```

Predict matches using the `predict` function. The output has a similar
structure to that of the `mec` function.

``` r
result_sup <- predict(model, df_1_new, df_2_new)
result_sup
#> The algorithm predicted 5 matches.
#> The first 5 predicted matches are:
#>        a     b ratio / 1000
#>    <num> <num>        <num>
#> 1:     3     3    0.6466869
#> 2:     4     4    0.5865049
#> 3:     1     1    0.5696382
#> 4:     5     5    0.3103742
#> 5:     2     2    0.2935612
#> ========================================================
#> The construction of the classification set was based on estimates of its size.
#> Estimated false link rate (FLR): 1.1486 %.
#> Estimated missing match rate (MMR): 1.1486 %.
```

## Funding

Work on this package is supported by the National Science Centre, OPUS
20 grant no. 2020/39/B/HS4/00941 (Towards census-like statistics for
foreign-born populations – quality, data integration and estimation).

## References

Lee, D., Zhang, L.-C. and Kim, J. K. (2022). [Maximum entropy
classification for record
linkage.](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm)
Survey Methodology, Statistics Canada, Catalogue No. 12-001-X, Vol. 48,
No. 1.

Vo, T. H., Chauvet, G., Happe, A., Oger, E., Paquelet, S., and Garès, V.
(2023). [Extending the Fellegi-Sunter record linkage model for
mixed-type data with application to the French national health data
system.](https://ideas.repec.org/a/eee/csdana/v179y2023ics0167947322002365.html)
Computational Statistics & Data Analysis, 179, 107656.

Sugiyama, M., Suzuki, T., Nakajima, S. et al. [Direct importance
estimation for covariate shift
adaptation.](https://doi.org/10.1007/s10463-008-0197-x) Ann Inst Stat
Math 60, 699–746 (2008).
