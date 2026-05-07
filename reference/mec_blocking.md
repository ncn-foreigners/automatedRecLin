# Blocked Unsupervised Maximum Entropy Classifier for Record Linkage

Runs graph-based blocking using
[blocking()](https://ncn-foreigners.ue.poznan.pl/blocking/reference/blocking.html),
fits one pooled unsupervised maximum entropy classifier (MEC) on
selected within-block pairs, and applies the fitted density-ratio model
blockwise.

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
  min_training_pairs = NULL,
  min_training_nonmatches = NULL,
  block_sampling_seed = NULL,
  nonmatch_sample_size = NULL,
  nonmatch_sampling_seed = NULL,
  prob_ratio = "2",
  start_params = NULL,
  nonpar_hurdle = TRUE,
  fixed_method = "Newton",
  delta = 0.5,
  eps = 0.05,
  max_iter_em = 10,
  tol_em = 1,
  controls_nleqslv = list(),
  controls_kliep = control_kliep(),
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

- min_training_pairs:

  Minimum number of within-block training pairs. If `NULL`,
  `min_training_nonmatches` should also be `NULL`, and all blocks are
  used for training.

- min_training_nonmatches:

  Minimum lower bound on within-block nonmatches. If `NULL`,
  `min_training_pairs` should also be `NULL`, and all blocks are used
  for training.

- block_sampling_seed:

  Optional seed for random training-block sampling.

- nonmatch_sample_size:

  Number of pairs sampled from the full Cartesian product of `A` and `B`
  to estimate nonmatch distribution parameters. If `NULL`, all pairs are
  used.

- nonmatch_sampling_seed:

  Optional seed for nonmatch pair sampling.

- prob_ratio:

  Probability/density ratio type (`"1"` or `"2"`). The default `"2"`
  uses the blockwise fixed-point equation.

- start_params:

  Start parameters for the `"binary"` and `"continuous_parametric"`
  methods.

- nonpar_hurdle:

  Currently unused in `mec_blocking()`.

- fixed_method:

  A method for solving blockwise fixed-point equations using the
  [FixedPoint()](https://rdrr.io/pkg/FixedPoint/man/FixedPoint.html)
  function.

- delta:

  A numeric value specifying the tolerance for the change in the
  estimated number of matches between MEC iterations.

- eps:

  A numeric value specifying the tolerance for the change in model
  parameters between MEC iterations.

- max_iter_em:

  Currently unused in `mec_blocking()`.

- tol_em:

  Currently unused in `mec_blocking()`.

- controls_nleqslv:

  Controls passed to the
  [nleqslv()](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html)
  function (only if the `"continuous_parametric"` method has been chosen
  for at least one variable).

- controls_kliep:

  Currently unused in `mec_blocking()`.

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

- `flr_est` – estimated false link rate (FLR),

- `mmr_est` – estimated missing match rate (MMR),

- `training_rule` – training-block selection rule used by the function,

- `block_estimates` – a `data.table` with block-level match-count and
  error-rate estimates,

- `training_blocks` – a `data.table` with blocks selected for pooled MEC
  training,

- `block_summary` – a `data.table` describing the final disjoint blocks,

- `excluded_records` – a list with records from `A` and `B` excluded by
  blocking,

- `pooled_model` – fitted pooled MEC density-ratio model used for
  blockwise scoring,

- `b_vars` – variables used for the `"binary"` method, with the prefix
  `"gamma_"`,

- `cpar_vars` – variables used for the `"continuous_parametric"` method,
  with the prefix `"gamma_"`,

- `cnonpar_vars` – variables used for the `"continuous_nonparametric"`
  method, currently `NULL`,

- `hm_vars` – variables used for the `"hit_miss"` method, currently
  `NULL`,

- `b_params` – parameters estimated using the `"binary"` method,

- `cpar_params` – parameters estimated using the
  `"continuous_parametric"` method,

- `cnonpar_params` – parameters estimated using the
  `"continuous_nonparametric"` method, currently `NULL`,

- `hm_params` – parameters estimated using the `"hit_miss"` method,
  currently `NULL`,

- `ratio_kliep` – result of
  [kliep()](https://thomvolker.github.io/densityratio/reference/kliep.html),
  currently `NULL`,

- `ratio_kliep_list` – variable-specific KLIEP results, currently
  `NULL`,

- `variables` – key variables used for comparison,

- `comparators` – comparison functions used to create comparison
  vectors,

- `methods` – MEC estimation methods used for the key variables,

- `nonmatch_sample_size` – number of full Cartesian-product pairs used
  to estimate nonmatch parameters,

- `nonmatch_sampling_seed` – seed used for nonmatch-pair sampling,

- `prob_ratio` – probability/density ratio type used for blockwise
  match-count estimation,

- `delta` – tolerance for changes in the estimated number of matches,

- `eps` – tolerance for changes in model parameters,

- `controls_nleqslv` – controls passed to
  [nleqslv()](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html),

- `controls_blocking` – additional arguments passed to
  [blocking()](https://ncn-foreigners.ue.poznan.pl/blocking/reference/blocking.html),

- `blocking_result` – raw object returned by
  [blocking()](https://ncn-foreigners.ue.poznan.pl/blocking/reference/blocking.html)
  if `keep_blocking_result = TRUE`; otherwise `NULL`,

- `training_Omega` – pooled training comparison vectors if
  `keep_training_data = TRUE`; otherwise `NULL`,

- `blocking_eval` – blocking diagnostics if `true_matches` is provided;
  otherwise `NULL`,

- `eval_metrics` – standard linkage quality metrics if `true_matches` is
  provided; otherwise `NULL`,

- `confusion` – confusion matrix if `true_matches` is provided;
  otherwise `NULL`.

## Details

The function assumes one-to-one linkage. The blocking stage defines
disjoint bipartite blocks. MEC is trained once on the pooled union of
selected within-block Cartesian products and is then applied separately
in each final block. The ANN distance returned by
[blocking()](https://ncn-foreigners.ue.poznan.pl/blocking/reference/blocking.html)
is not used.

If both `min_training_pairs` and `min_training_nonmatches` are `NULL`,
all final blocks are used for pooled training. If both are supplied,
blocks are sampled without replacement until both thresholds are
reached. Supplying only one threshold is an error.

Nonmatch distribution parameters are estimated from a simple random
sample from the full Cartesian product of `A` and `B`, not from the
selected training blocks.

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
#> Blocked MEC record linkage based on the following variables:  
#> name, surname, city.
#> ========================================================
#> Number of final blocks: 5.
#> Training rule: all_blocks.
#> Number of training blocks: 5.
#> Number of training pairs: 5.
#> Training nonmatch lower bound: 0.
#> ========================================================
#> The algorithm predicted 5 matches.
#> The first 5 predicted matches are:
#>        a     b block ratio / 1000
#>    <int> <int> <num>        <num>
#> 1:     1     1     1   0.06944444
#> 2:     2     2     2   0.06944444
#> 3:     3     3     3   0.06944444
#> 4:     4     4     4   0.06944444
#> 5:     5     5     5   0.06944444
#> ========================================================
#> Estimated false link rate (FLR): 0.0000 %.
#> Estimated missing match rate (MMR): 0.0000 %.
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
