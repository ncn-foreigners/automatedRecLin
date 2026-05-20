# Blocked Unsupervised Maximum Entropy Classifier for Record Linkage

Runs graph-based blocking using
[blocking()](https://ncn-foreigners.ue.poznan.pl/blocking/reference/blocking.html),
defines a blocking candidate-pair space, and fits an inverted
unsupervised maximum entropy classifier (MEC) directly on all candidate
pairs.

## Usage

``` r
mec_blocking(
  A,
  B,
  variables,
  comparators = NULL,
  methods = NULL,
  blocking_x = NULL,
  blocking_y = NULL,
  blocking_variables = variables,
  blocking_sep = " ",
  controls_blocking = list(),
  start_params = NULL,
  alpha = 0,
  delta = 0.5,
  eps = 0.05,
  controls_nleqslv = list(),
  true_matches = NULL,
  keep_blocking_result = FALSE,
  keep_training_data = FALSE,
  verbose = FALSE
)
```

## Arguments

- A:

  A duplicate-free `data.frame` or `data.table`.

- B:

  A duplicate-free `data.frame` or `data.table`.

- variables:

  A character vector of key variables used to create MEC comparison
  vectors.

- comparators:

  A named list of functions for comparing pairs of records.

- methods:

  A named list of methods used for estimation (`"binary"` or
  `"continuous_parametric"`). Other unsupervised MEC methods are not
  supported by `mec_blocking()` at this stage.

- blocking_x:

  Optional input passed as `x` to
  [blocking()](https://ncn-foreigners.ue.poznan.pl/blocking/reference/blocking.html).

- blocking_y:

  Optional input passed as `y` to
  [blocking()](https://ncn-foreigners.ue.poznan.pl/blocking/reference/blocking.html).

- blocking_variables:

  Variables used to create blocking strings when `blocking_x` and
  `blocking_y` are not supplied.

- blocking_sep:

  Separator used when concatenating `blocking_variables`.

- controls_blocking:

  A list of additional arguments passed to
  [blocking()](https://ncn-foreigners.ue.poznan.pl/blocking/reference/blocking.html),
  except `x` and `y`.

- start_params:

  Start parameters for the `"binary"` and `"continuous_parametric"`
  methods.

- alpha:

  A single numeric value in `[0, 1)` controlling the fraction of the
  current nonmatch complement dropped from nonmatch-side parameter
  fitting after the first inverted MEC iteration. The first U-side fit
  uses the full initial fitting set, and posterior/count formulas
  continue to use the full current nonmatch count.

- delta:

  A numeric value specifying the tolerance for the change in the
  estimated number of nonmatches between MEC iterations.

- eps:

  A numeric value specifying the tolerance for the change in model
  parameters between MEC iterations.

- controls_nleqslv:

  Controls passed to the
  [nleqslv()](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html)
  function (only if the `"continuous_parametric"` method has been chosen
  for at least one variable).

- true_matches:

  A `data.frame` or `data.table` indicating known matches.

- keep_blocking_result:

  Logical indicating whether to store the raw object returned by
  [blocking()](https://ncn-foreigners.ue.poznan.pl/blocking/reference/blocking.html).

- keep_training_data:

  Logical indicating whether to store pooled training comparison
  vectors.

- verbose:

  Logical indicating whether to print progress messages.

## Value

Returns a list of class `"mec_blocking"` containing:

- `M_est` – a `data.table` with predicted matches and columns `a`, `b`,
  `block`, and `ratio`,

- `n_M_est` – estimated total number of matches across all blocks,

- `n_U_est` – estimated total number of candidate nonmatches,

- `alpha` – fraction of the current nonmatch complement dropped from
  later U-side fitting,

- `candidate_pair_count` – number of candidate pairs in \\\Omega_B\\,

- `block_estimates` – a `data.table` with block-level size and
  match-count diagnostics,

- `block_summary` – a `data.table` describing the final disjoint blocks,

- `excluded_records` – a list with records from `A` and `B` excluded by
  blocking,

- `b_vars` – variables used for the `"binary"` method, with the prefix
  `"gamma_"`,

- `cpar_vars` – variables used for the `"continuous_parametric"` method,
  with the prefix `"gamma_"`,

- `b_params` – parameters estimated using the `"binary"` method,

- `cpar_params` – parameters estimated using the
  `"continuous_parametric"` method,

- `variables` – key variables used for comparison,

- `comparators` – comparison functions used to create comparison
  vectors,

- `methods` – MEC estimation methods used for the key variables,

- `delta` – tolerance for changes in the estimated number of nonmatches,

- `eps` – tolerance for changes in nonmatch-side model parameters,

- `controls_nleqslv` – controls passed to
  [nleqslv()](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html),

- `blocking_result` – raw object returned by
  [blocking()](https://ncn-foreigners.ue.poznan.pl/blocking/reference/blocking.html)
  if `keep_blocking_result = TRUE`; otherwise `NULL`,

- `training_Omega` – candidate-space comparison vectors with inverted
  scores if `keep_training_data = TRUE`; otherwise `NULL`,

- `blocking_eval` – blocking diagnostics if `true_matches` is provided;
  otherwise `NULL`,

- `mec_eval` – MEC-selection diagnostics among known matches retained in
  the candidate-pair space if `true_matches` is provided; otherwise
  `NULL`,

- `eval_metrics` – empirical linkage quality metrics based on
  `true_matches`; otherwise `NULL`,

- `confusion` – empirical confusion matrix based on `true_matches`;
  otherwise `NULL`.

## Details

The function assumes one-to-one linkage. The blocking stage defines
disjoint bipartite blocks, and the candidate-pair space \\\Omega_B\\ is
the union of within-block Cartesian products. Duplicate candidate pairs
are removed deterministically before MEC fitting.

The blocked MEC fit is inverted relative to
[`mec()`](https://ncn-foreigners.github.io/automatedRecLin/reference/mec.md).
The initial match set contains at most \\\nu\\ feasible pairs, where
\\\nu\\ is the structural one-to-one upper bound. Initial feasible
matches are selected greedily by an unweighted disagreement norm: binary
agreement indicators use `1 - gamma`, while continuous dissimilarities
use `gamma` unchanged. At each iteration, match-side parameters are
estimated from the current greedy one-to-one match set, and
nonmatch-side parameters are estimated from its complement.

The `alpha` argument applies only to nonmatch-side distribution
estimation. The first U-side fit uses the full initial complement. In
later iterations, the least reliable current nonmatches are dropped from
the U-side fitting sample, with reliability ranked by the previous
nonmatch posterior estimate and then by the inverted density ratio if
the posterior is unavailable. The posterior and count updates still use
the full current complement size, and the final match set remains
one-to-one.

The returned `ratio` is \\s = u / m\\, where \\u\\ and \\m\\ denote the
estimated nonmatch and match comparison-vector densities. Smaller values
are therefore more match-like. Updated match sets are selected greedily
in ascending order of this ratio.

If \\N = \|\Omega_B\|\\ and \\\nu\\ is the maximum feasible one-to-one
matching size in the candidate graph, the estimated number of nonmatches
is bounded below by \\N - \nu\\. For the disjoint complete blocks
reconstructed by this function, \\\nu = \sum_h \min(n\_{Ah}, n\_{Bh})\\.

If the initialized match set exhausts the candidate-pair space, for
example when \\N = \nu\\, there is no candidate complement from which to
estimate nonmatch parameters. In that case the function returns the
structurally feasible initialized match set, sets `n_U_est = 0`, and
leaves nonmatch-side parameters unavailable.

## Author

Adam Struzik

## Examples

``` r
df_1 <- data.frame(
  name = c("Emma", "Liam", "Olivia", "Noah", "Ava"),
  surname = c("Smith", "Jones", "Brown", "Davis", "Miller"),
  city = c("Boston", "Boston", "Austin", "Austin", "Denver")
)
df_2 <- data.frame(
  name = c("Emma", "Liam", "Olivia", "Noah", "Ava"),
  surname = c("Smith", "Jones", "Brown", "Davis", "Miller"),
  city = c("Boston", "Boston", "Austin", "Austin", "Denver")
)

blocking_x <- matrix(
  c(1, 0, 0, 1, 1, 1, 2, 0, 0, 2),
  ncol = 2,
  byrow = TRUE
)
blocking_y <- blocking_x

result <- mec_blocking(
  A = df_1,
  B = df_2,
  variables = c("name", "surname", "city"),
  blocking_x = blocking_x,
  blocking_y = blocking_y,
  controls_blocking = list(
    representation = "custom_matrix",
    ann = "kd",
    distance = "euclidean",
    seed = 1
  ),
  true_matches = data.frame(a = 1:5, b = 1:5)
)
result
#> Blocked MEC record linkage based on:  
#> name, surname, city.
#> ========================================================
#> The algorithm predicted 5 matches.
#> The first 5 predicted matches are:
#>        a     b block ratio / 1000
#>    <int> <int> <num>        <num>
#> 1:     1     1     1            0
#> 2:     2     2     2            0
#> 3:     3     3     3            0
#> 4:     4     4     4            0
#> 5:     5     5     5            0
#> ========================================================
#> ========================================================
#> Blocking diagnostics:
#> Known matches: 5.
#> Known matches retained by blocking: 5.
#> Known matches missed by blocking: 0.
#> Blocking MMR: 0.0000 %.
#> Candidate pairs retained: 5 of 25.
#> Candidate pair reduction: 80.0000 %.
#> ========================================================
#> Evaluation metrics:
#> FLR (%) MMR (%) 
#>  0.0000  0.0000 
```
