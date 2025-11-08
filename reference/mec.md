# Unsupervised Maximum Entropy Classifier for Record Linkage

Implements several extensions to the maximum entropy classification
(MEC) algorithm for record linkage (see [Lee et al.
(2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm)),
iteratively estimating probability/density ratios to classify record
pairs into matches and non-matches based on comparison vectors.

## Usage

``` r
mec(
  A,
  B,
  variables,
  comparators = NULL,
  methods = NULL,
  duplicates_in_A = FALSE,
  start_params = NULL,
  nonpar_hurdle = TRUE,
  set_construction = NULL,
  target_rate = 0.03,
  max_iter_bisection = 100,
  tol = 0.005,
  delta = 0.5,
  eps = 0.05,
  max_iter_em = 10,
  tol_em = 1,
  controls_nleqslv = list(),
  controls_kliep = control_kliep(),
  true_matches = NULL
)
```

## Arguments

- A:

  A duplicate-free `data.frame` or `data.table`.

- B:

  A duplicate-free `data.frame` or `data.table`.

- variables:

  A character vector of key variables used to create comparison vectors.

- comparators:

  A named list of functions for comparing pairs of records.

- methods:

  A named list of methods used for estimation (`"binary"`,
  `"continuous_parametric"`, `"continuous_nonparametric"` or
  `"hit_miss"`).

- duplicates_in_A:

  Logical indicating whether to allow `A` to have duplicate records.

- start_params:

  Start parameters for the `"binary"`, `"continuous_parametric"` and
  `"hit_miss"` methods.

- nonpar_hurdle:

  Logical indicating whether to use a hurdle model or not (used only if
  the `"continuous_nonparametric"` method has been chosen for at least
  one variable).

- set_construction:

  A method for constructing the predicted set of matches (`"size"`,
  `"flr"` or `"mmr"`).

- target_rate:

  A target false link rate (FLR) or missing match rate (MMR) (used only
  if `set_construction == "flr"` or `set_construction == "mmr"`).

- max_iter_bisection:

  A maximum number of iterations for the bisection procedure (used only
  if `set_construction == "flr"` or `set_construction == "mmr"`).

- tol:

  Error tolerance in the bisection procedure (used only if
  `set_construction == "flr"` or `set_construction == "mmr"`).

- delta:

  A numeric value specifying the tolerance for the change in the
  estimated number of matches between iterations.

- eps:

  A numeric value specifying the tolerance for the change in model
  parameters between iterations.

- max_iter_em:

  A maximum number of iterations for the EM algorithm (used only if the
  `"hit_miss"` method has been chosen for at least one variable).

- tol_em:

  Error tolerance in the EM algorithm (used only if the `"hit_miss"`
  method has been chosen for at least one variable).

- controls_nleqslv:

  Controls passed to the
  [nleqslv](https://rdrr.io/pkg/nleqslv/man/nleqslv.html) function (only
  if the `"continuous_parametric"` method has been chosen for at least
  one variable).

- controls_kliep:

  Controls passed to the
  [kliep](https://thomvolker.github.io/densityratio/reference/kliep.html)
  function (only if the `"continuous_nonparametric"` method has been
  chosen for at least one variable).

- true_matches:

  A `data.frame` or `data.table` indicating known matches.

## Value

Returns a list containing:  

- `M_est` – a `data.table` with predicted matches,

- `n_M_est` – estimated classification set size,

- `flr_est` – estimated false link rate (FLR),

- `mmr_est` – estimated missing match rate (MMR),

- `iter_bisection` – the number of iterations in the bisection
  procedure,

- `b_vars` – a character vector of variables used for the `"binary"`
  method (with the prefix `"gamma_"`),

- `cpar_vars` – a character vector of variables used for the
  `"continuous_parametric"` method (with the prefix `"gamma_"`),

- `cnonpar_vars` – a character vector of variables used for the
  `"continuous_nonparametric"` method (with the prefix `"gamma_"`),

- `hm_vars` – a character vector of variables used for the `"hit_miss"`
  method (with the prefix `"gamma_"`),

- `b_params` – parameters estimated using the `"binary"` method,

- `cpar_params` – parameters estimated using the
  `"continuous_parametric"` method,

- `hm_params` – parameters estimated using the `"hit_miss"` method,

- `ratio_kliep` – a result of the
  [kliep](https://thomvolker.github.io/densityratio/reference/kliep.html)
  function,

- `variables` – a character vector of key variables used for comparison,

- `set_construction` – a method for constructing the predicted set of
  matches,

- `eval_metrics` – standard metrics for quality assessment (if
  `true_matches` is provided),

- `confusion` – confusion matrix (if `true_matches` is provided).

## Details

Consider two datasets without duplicates: \\A\\ and \\B\\. Let the
bipartite comparison space \\\Omega = A \times B\\ consist of matches
\\M\\ and non-matches \\U\\ between the records in files \\A\\ and
\\B\\. For any pair of records \\(a,b) \in \Omega\\, let
\\\pmb{\gamma}\_{ab} = (\gamma\_{ab}^1,\gamma\_{ab}^2,
\ldots,\gamma\_{ab}^K)'\\ be the comparison vector between a set of key
variables. The original MEC algorithm uses the binary comparison
function to evaluate record pairs across two datasets. However, this
approach may be insufficient when handling datasets with frequent errors
across multiple variables.

We propose the use of continuous comparison functions to address the
limitations of binary comparison methods. We consider every semi-metric,
i.e., a function \\d: A \times B \to \mathbb{R}\\, satisfying the
following conditions:  

1.  \\d(x,y) \geq 0\\,

2.  \\d(x,y) = 0\\ if and only if \\x = y\\,

3.  \\d(x,y) = d(y,x)\\.

For example, we can use \\1 - \text{Jaro-Winkler distance}\\ for
character variables (which is implemented in the `automatedRecLin`
package as the `jarowinkler_complement` function) or the Euclidean
distance for numerical variables. The `automatedRecLin` package allows
the use of a different comparison function for each key variable (which
should be specified as a list in the `comparators` argument). The
default function for each key variable is
[cmp_identical](https://rdrr.io/pkg/reclin2/man/comparators.html) (the
binary comparison function).

The `mec` function offers different approaches to estimate the
probability/density ratio between matches and non-matches, which should
be specified as a list in the `methods` argument. The available methods
suitable for the binary comparison function are `"binary"` and
`"hit_miss"`. Both assume that \\\gamma\_{ab}^k\|M\\ and
\\\gamma\_{ab}^k\|U\\ follow Bernoulli distributions. `"binary"` and
`"hit_miss"` both estimate the parameters for the matches iteratively,
but `"binary"` estimates the parameters for the non-matches only at the
start, while `"hit_miss"` does so iteratively using a hit-miss model
(for details see [Lee et al.
(2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm)).
`"binary"` is the default method for each variable.

For the continuous semi-metrics we suggest the usage of
`"continuous_parametric"` or `"continuous_nonparametric"` method. The
`"continuous_parametric"` method assumes that \\\gamma\_{ab}^k\|M\\ and
\\\gamma\_{ab}^k\|U\\ follow hurdle Gamma distributions. The density
function of a hurdle Gamma distribution is characterized by three
parameters \\p_0 \in (0,1)\\ and \\\alpha, \beta \> 0\\ as follows: \$\$
f(x;p_0,\alpha,\beta) = p_0^{\mathbb{I}(x = 0)}\[(1 - p_0)
v(x;\alpha,\beta)\]^{\mathbb{I}(x \> 0)}, \$\$ where \$\$
v(x;\alpha,\beta) = \frac{\beta^{\alpha} x^{\alpha - 1} \exp(-\beta x)}
{\Gamma(\alpha)} \$\$ is the density function of a Gamma distribution
(for details see [Vo et al.
(2023)](https://ideas.repec.org/a/eee/csdana/v179y2023ics0167947322002365.html)).
At the beginning, the algorithm estimates the parameters for the
non-matches and then does it iteratively for the matches. The
`"continuous_nonparametric"` method does not assume anything about the
distributions of the comparison vectors. It iteratively directly
estimates the density ratio between the matches and the non-matches,
using the Kullback-Leibler Importance Estimation Procedure (KLIEP). For
details see [Sugiyama et al.
(2008)](https://link.springer.com/article/10.1007/s10463-008-0197-x).

The `mec` function allows the construction of the predicted set of
matches using its estimated size or the bisection procedure, described
in [Lee et al.
(2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm),
based on a target False Link Rate (FLR) or missing match rate (MMR). To
use the second option, set `set_construction = "flr"` or
`set_construction = "mmr"` and specify a target error rate using the
`target_rate` argument.

The assumption that \\A\\ and \\B\\ contain no duplicate records might
be relaxed by allowing \\A\\ to have duplicates. To do so, set
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
true_matches <- data.frame(
  "a" = 1:8,
  "b" = 1:8
)

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

set.seed(1)
result <- mec(A = df_1, B = df_2,
              variables = variables,
              comparators = comparators,
              methods = methods,
              true_matches = true_matches)
result
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
#> ========================================================
#> Evaluation metrics:
#> FLR MMR 
#>   0   0 
```
