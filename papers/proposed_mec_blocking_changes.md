# Proposed changes to `mec_blocking()`

This document describes two implementation changes for `mec_blocking()`.

## 1. Classify singleton blocks using only the density ratio

A singleton block is a block `c` with exactly one record from file `A` and exactly one record from file `B`:

```text
n_A_c = 1
n_B_c = 1
```

Hence the block contains exactly one candidate pair:

```text
Omega_c = { (a, b) }
```

### Required behavior

For singleton blocks, do not use the blockwise fixed-point equation and do not estimate the local number of matches through the usual MEC count-estimation procedure.

Instead, classify the only candidate pair directly using its fitted density ratio:

```text
r_ab_c = f_M(gamma_ab) / f_U(gamma_ab)
```

Use the rule:

```text
if r_ab_c >= 1:
    classify (a, b) as a match
else:
    classify (a, b) as a nonmatch
```

Equivalently:

```text
n_M_c_hat = 1 if r_ab_c >= 1
n_M_c_hat = 0 if r_ab_c < 1
```

and

```text
M_c_hat = { (a, b) } if r_ab_c >= 1
M_c_hat = empty set otherwise
```

### Rationale

For a singleton block, there is no competition between candidate pairs and no one-to-one conflict to resolve. The task is a binary classification problem for one pair.

The rule `r_ab_c >= 1` means that the observed comparison vector is at least as likely under the fitted match distribution as under the fitted nonmatch distribution.

This also avoids the degenerate behavior of the current fixed-point approach, where `1` is always a fixed-point solution for singleton blocks.

### Edge cases

Handle edge cases as follows:

```text
r_ab_c >= 1       -> match
0 <= r_ab_c < 1   -> nonmatch
r_ab_c = Inf      -> match
r_ab_c = NA/NaN   -> nonmatch
```

If the implementation works on the log scale, use:

```text
log(r_ab_c) >= 0
```

For all non-singleton blocks, keep the existing blockwise MEC procedure.

## 2. Remove model-based error-rate estimators

The function should no longer compute or return model-based estimates of FLR and MMR.

Remove or stop returning the following global output components:

```text
flr_est
mmr_est
```

Also remove block-level estimated error-rate columns from `block_estimates`.

### Rationale

The model-based FLR/MMR estimates are not sufficiently reliable in the blocked MEC setting and should not be presented as standard output.

The main purpose of `mec_blocking()` should be to return predicted links, fitted density ratios, and basic block diagnostics, not model-based error-rate estimates.

This also removes the need to assign posterior-like matching probabilities to singleton pairs.

### Output after the change

The returned object should focus on linkage results and diagnostics such as:

```text
M_est
n_M_est
training_rule
training_blocks
block_estimates
block_summary
excluded_records
pooled_model
blocking_eval
eval_metrics
confusion
```

`M_est` should still contain predicted links and their fitted density ratios:

```text
a
b
block
ratio
```

`block_estimates` should keep block-size and match-count information, for example:

```text
block
n_A_c
n_B_c
N_c
L_c
n_M_c_hat
n_selected_c
```

but it should not contain estimated FLR/MMR columns.

### Known matches

If `known_matches` is supplied, empirical evaluation metrics may still be reported.

These are not model-based error-rate estimates. They are direct evaluation metrics computed by comparing predicted links with the supplied known matches.

Therefore, the following optional outputs may remain when `known_matches` is provided:

```text
blocking_eval
eval_metrics
confusion
```

## Summary

For each block:

```text
if n_A_c == 1 and n_B_c == 1:
    compute r_ab_c
    if r_ab_c >= 1:
        select the pair
    else:
        select no pair
else:
    use the existing blockwise MEC procedure
```

The function should no longer compute or return model-based FLR/MMR estimates.

Empirical evaluation against `known_matches`, if supplied, should remain available.
