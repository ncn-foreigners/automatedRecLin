# MEC with blocking

## Setup

Load required packages.

``` r

library(automatedRecLin)
library(data.table)

options("text2vec.mc.cores" = 1L)
```

## Data

We use the full example Census and Customer Information System (CIS)
datasets from [McLeod et
al. (2011)](https://wayback.archive-it.org/12090/20231221144450/https://cros-legacy.ec.europa.eu/content/job-training_en).
The goal is to link records from CIS to records from Census.

``` r

data("census", package = "automatedRecLin")
data("cis", package = "automatedRecLin")
setDT(census)
setDT(cis)

NROW(cis)
#> [1] 24613
NROW(census)
#> [1] 25343
```

The `person_id` variable identifies the correct linkage. We use this
information only to evaluate the result.

``` r

cis[is.na(cis)] <- ""
census[is.na(census)] <- ""

cis[, pername1 := gsub("-", "", pername1)]
census[, pername1 := gsub("-", "", pername1)]

true_matches <- merge(
  x = cis[, .(a = .I, person_id)],
  y = census[, .(b = .I, person_id)],
  by = "person_id"
)[, .(a, b)]

NROW(true_matches)
#> [1] 24043
```

## MEC with blocking

We compare forename and surname using the Jaro-Winkler distance. These
two comparison variables are modeled with the continuous parametric MEC
method. Sex and date-of-birth variables use the default binary method.
Address fields are used only to construct blocks.

``` r

variables <- c(
  "pername1", "pername2", "sex",
  "dob_day", "dob_mon", "dob_year"
)

comparators <- list(
  "pername1" = jarowinkler_complement(),
  "pername2" = jarowinkler_complement()
)

methods <- list(
  "pername1" = "continuous_parametric",
  "pername2" = "continuous_parametric"
)

blocking_variables <- c(variables, "enumcap", "enumpc")
```

Run blocked MEC. The model is trained on all candidate pairs retained by
blocking.

``` r

set.seed(1)

result <- mec_blocking(
  A = cis,
  B = census,
  variables = variables,
  comparators = comparators,
  methods = methods,
  blocking_variables = blocking_variables,
  blocking_sep = " ",
  controls_blocking = list(seed = 1, n_threads = 1),
  alpha = 0.5,
  true_matches = true_matches
)

result
#> Blocked MEC record linkage based on:  
#> pername1, pername2, sex, dob_day, dob_mon, dob_year.
#> ========================================================
#> The algorithm predicted 23697 matches.
#> The first 6 predicted matches are:
#>        a     b block ratio / 1000
#>    <int> <int> <num>        <num>
#> 1: 12264 18361 17172 3.428838e-12
#> 2: 23367 13031 12223 3.428838e-12
#> 3: 23495 15194 14243 3.428838e-12
#> 4:  1768 12279 11529 7.032757e-12
#> 5:   343 18657 17446 8.007616e-12
#> 6:  2124  5497  5152 8.007616e-12
#> ========================================================
#> ========================================================
#> Blocking diagnostics:
#> Known matches: 24043.
#> Known matches retained by blocking: 23687.
#> Known matches missed by blocking: 356.
#> Blocking MMR: 1.4807 %.
#> Candidate pairs retained: 25343 of 623767259.
#> Candidate pair reduction: 99.9959 %.
#> ========================================================
#> Evaluation metrics:
#> FLR (%) MMR (%) 
#>  0.0506  1.4890
```

## Blocking efficiency and linkage results

The full Cartesian product contains 623,767,259 record pairs. Blocking
reduces this to 25,343 candidate pairs, while retaining 98.52% of known
links. The final linkage set contains 23,697 predicted matches.

    #>        step                                result
    #>      <char>                                <char>
    #> 1: Training  all_candidate_pairs on 23,725 blocks
    #> 2: Blocking 23,687 of 24,043 known links retained
    #> 3:  Linkage              FLR = 0.05%; MMR = 1.49%
