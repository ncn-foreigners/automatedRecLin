# Train a Record Linkage Model

Trains a supervised record linkage model using probability or density
ratio estimation, based on [Lee et al.
(2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm),
with several extensions.

## Usage

``` r
train_rec_lin(
  A,
  B,
  matches,
  variables,
  comparators = NULL,
  methods = NULL,
  prob_ratio = NULL,
  nonpar_hurdle = TRUE,
  controls_nleqslv = list(),
  controls_kliep = control_kliep()
)
```

## Arguments

- A:

  A duplicate-free `data.frame` or `data.table`.

- B:

  A duplicate-free `data.frame` or `data.table`.

- matches:

  A `data.frame` or `data.table` indicating known matches.

- variables:

  A character vector of key variables used to create comparison vectors.

- comparators:

  A named list of functions for comparing pairs of records.

- methods:

  A named list of methods used for estimation (`"binary"`,
  `"continuous_parametric"` or `"continuous_nonparametric"`).

- prob_ratio:

  Probability/density ratio type (`"1"` or `"2"`).

- nonpar_hurdle:

  Logical indicating whether to use a hurdle model or not (used only if
  the `"continuous_nonparametric"` method has been chosen for at least
  one variable).

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

## Value

Returns a list containing:  

- `b_vars` – a character vector of variables used for the `"binary"`
  method (with the prefix `"gamma_"`),

- `cpar_vars` – a character vector of variables used for the
  `"continuous_parametric"` method (with the prefix `"gamma_"`),

- `cnonpar_vars` – a character vector of variables used for the
  `"continuous_nonparametric"` method (with the prefix `"gamma_"`),

- `b_params` – parameters estimated using the `"binary"` method,

- `cpar_params` – parameters estimated using the
  `"continuous_parametric"` method,

- `cnonpar_params` – probability of exact matching estimated using the
  `"continuous_nonparametric"` method,

- `ratio_kliep` – a result of the
  [kliep](https://thomvolker.github.io/densityratio/reference/kliep.html)
  function,

- `ratio_kliep_list` – an object containing the results of the
  [kliep](https://thomvolker.github.io/densityratio/reference/kliep.html)
  function,

- `ml_model` – here `NULL`,

- `pi_est` – a prior probability of matching,

- `match_prop` – proportion of matches in the smaller dataset,

- `variables` – a character vector of key variables used for comparison,

- `comparators` – a list of functions used to compare pairs of records,

- `methods` – a list of methods used for estimation,

- `"prob_ratio"` – probability/density ratio type.

## Details

Consider two datasets: \\A\\ and \\B\\. Let the bipartite comparison
space \\\Omega = A \times B\\ consist of matches \\M\\ and non-matches
\\U\\ between the records in files \\A\\ and \\B\\. For any pair of
records \\(a,b) \in \Omega\\, let \\\pmb{\gamma}\_{ab} =
(\gamma\_{ab}^1,\gamma\_{ab}^2, \ldots,\gamma\_{ab}^K)'\\ be the
comparison vector between a set of key variables. The original MEC
algorithm uses the binary comparison function to evaluate record pairs
across two datasets. However, this approach may be insufficient when
handling datasets with frequent errors across multiple variables.

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

The `train_rec_lin` function is used to train a record linkage model,
when \\M\\ and \\U\\ are known (which might later serve as a classifier
for pairs outside \\\Omega\\). It offers different approaches to
estimate the probability/density ratio between matches and non-matches,
which should be specified as a list in the methods argument. The method
suitable for the binary comparison function is `"binary"`, which is also
the default method for each variable.

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
The `"continuous_nonparametric"` method does not assume anything about
the distributions of the comparison vectors. It directly estimates the
density ratio between the matches and the non-matches, using the
Kullback-Leibler Importance Estimation Procedure (KLIEP). For details
see [Sugiyama et al.
(2008)](https://link.springer.com/article/10.1007/s10463-008-0197-x).

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
model
#> Record linkage model based on the following variables: name, surname.
#> The prior probability of matching is 0.04.
#> ========================================================
#> Variables selected for the continuous nonparametric method: name, surname.
#> ========================================================
#> Probability/density ratio type: 1.
```
