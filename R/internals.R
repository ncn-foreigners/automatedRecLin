# Avoid CRAN notes for data.table column names used through non-standard evaluation.
utils::globalVariables(c(
  ".",
  "a",
  "b",
  "block",
  "cumulative_nonmatches_min",
  "cumulative_pairs",
  "n_A",
  "n_B",
  "nonmatches_min",
  "pair_count",
  "selected_pairs"
))

#' @noRd
binary_formula <- function(df) {
  if (NCOL(df) == 0L) {
    return(numeric())
  }

  colMeans(as.matrix(df))
}

#' @noRd
p_0_formula <- function(df) {
  if (NCOL(df) == 0L) {
    return(numeric())
  }

  df_matrix <- as.matrix(df)
  colSums(df_matrix == 0) / NROW(df_matrix)
}

#' @noRd
f_alpha <- function(alpha, gamma) {
  gamma <- gamma[gamma > 0]
  gamma_mean <- mean(gamma)
  sum(log(gamma) - log(gamma_mean) - digamma(alpha) + log(alpha))
}

#' @noRd
f_alpha_iterative <- function(alpha, beta, gamma) {
  gamma <- gamma[gamma > 0]
  sum(log(gamma) + log(beta) - digamma(alpha))
}

#' @noRd
gamma_plus_formula <- function(df) {
  if (NCOL(df) == 0L) {
    return(numeric())
  }

  vapply(df, function(x) {
    x <- x[x > 0]
    mean(x)
  }, numeric(1))
}

#' @noRd
alpha_formula <- function(df, fun) {
  if (NCOL(df) == 0L) {
    return(numeric())
  }

  vapply(df, function(y) {
    fun(x = 1, fn = f_alpha, gamma = y, method = "Newton")$x
  }, numeric(1))
}

#' @noRd
alpha_formula_summary <- function(n_positive, positive_sum, log_positive_sum, fun) {
  vapply(seq_along(n_positive), function(i) {
    if (n_positive[i] < 2L) {
      stop("The nonmatch sample should contain at least two positive continuous comparisons for each continuous parametric variable.")
    }

    positive_mean <- positive_sum[i] / n_positive[i]
    f_alpha_summary <- function(alpha) {
      log_positive_sum[i] - n_positive[i] * log(positive_mean) -
        n_positive[i] * digamma(alpha) + n_positive[i] * log(alpha)
    }

    fun(x = 1, fn = f_alpha_summary, method = "Newton")$x
  }, numeric(1))
}

#' @noRd
alpha_formula_iterative <- function(df, fun, beta) {
  K <- length(beta)
  par <- sapply(1:K, function(x) {
    fun(x = 1, fn = f_alpha_iterative, gamma = df[, as.numeric(x)], beta = beta[x], method = "Newton")$x
  })
  par
}

#' @importFrom stats dgamma
#'
#' @noRd
hurdle_gamma_density <- function(x, p_0, alpha, beta) {
  ifelse(x == 0, p_0, 1) *
    ifelse(x > 0, (1 - p_0) * stats::dgamma(x = x, shape = alpha, scale = 1 / beta), 1)
}

#' @noRd
# Check for any exact agreement using only key variables to reduce copying.
has_perfect_agreement <- function(A, B, variables) {
  A_keys <- unique(A[, variables, with = FALSE])
  B_keys <- unique(B[, variables, with = FALSE])

  NROW(A_keys[B_keys, on = variables, nomatch = 0L, mult = "first"]) > 0
}

#' @noRd
positive_gamma_table <- function(gamma, variable) {
  gamma_dt <- data.table::data.table(value = gamma[gamma > 0])
  data.table::setnames(gamma_dt, "value", variable)
  gamma_dt
}

#' @noRd
validate_match_pairs <- function(matches, n_A, n_B, arg_name = "matches") {
  if (!(is.data.frame(matches) || is.data.table(matches))) {
    stop(sprintf("`%s` should be a data.frame or a data.table.", arg_name))
  }

  data.table::setDT(matches)

  if (!(length(colnames(matches)) == 2L && all(colnames(matches) == c("a", "b")))) {
    stop(sprintf("`%s` should consist of two columns: a, b.", arg_name))
  }

  if (anyNA(matches[["a"]]) || anyNA(matches[["b"]])) {
    stop(sprintf("`%s` cannot contain missing values.", arg_name))
  }

  if (!is.numeric(matches[["a"]]) || !is.numeric(matches[["b"]])) {
    stop(sprintf("`%s` should contain numeric row indices.", arg_name))
  }

  if (any(matches[["a"]] < 1) || any(matches[["b"]] < 1)) {
    stop(sprintf("`%s` should contain positive row indices.", arg_name))
  }

  if (any(matches[["a"]] != as.integer(matches[["a"]])) ||
      any(matches[["b"]] != as.integer(matches[["b"]]))) {
    stop(sprintf("`%s` should contain integer row indices.", arg_name))
  }

  if (any(matches[["a"]] > n_A) || any(matches[["b"]] > n_B)) {
    stop(sprintf("`%s` contains row indices outside the input datasets.", arg_name))
  }

  if (anyDuplicated(matches[, c("a", "b"), with = FALSE]) > 0L) {
    stop(sprintf("`%s` should not contain duplicate record pairs.", arg_name))
  }

  invisible(matches)
}

#' @noRd
validate_choice <- function(value, choices, arg_name) {
  if (length(value) != 1L || is.na(value) || !(value %in% choices)) {
    stop(sprintf("`%s` should be one of: %s.", arg_name, paste(choices, collapse = ", ")))
  }

  value
}

#' @noRd
validate_methods <- function(methods, variables, allowed_methods, default_method = "binary") {
  if (is.null(methods)) {
    methods <- list()
  } else if (!is.list(methods)) {
    stop("`methods` should be a list.")
  }

  methods <- unlist(methods, use.names = TRUE)
  missing_variables <- variables[!(variables %in% names(methods))]
  methods[missing_variables] <- default_method
  methods <- methods[variables]

  invalid_methods <- unique(methods[!(methods %in% allowed_methods)])
  if (length(invalid_methods) > 0L) {
    stop(
      sprintf(
        "`methods` contains unsupported values: %s. Allowed methods are: %s.",
        paste(invalid_methods, collapse = ", "),
        paste(allowed_methods, collapse = ", ")
      )
    )
  }

  methods
}

#' @noRd
extract_method_variables <- function(methods, include_hit_miss = FALSE) {
  gamma_names_for <- function(method_name) {
    variable_names <- names(methods[methods == method_name])

    if (length(variable_names) == 0L) {
      character()
    } else {
      paste0("gamma_", variable_names)
    }
  }

  method_variables <- list(
    b_vars = gamma_names_for("binary"),
    cpar_vars = gamma_names_for("continuous_parametric"),
    cnonpar_vars = gamma_names_for("continuous_nonparametric")
  )

  if (include_hit_miss) {
    method_variables[["hm_vars"]] <- gamma_names_for("hit_miss")
  }

  method_variables
}

#' @noRd
drop_constant_key_variables <- function(A, B, variables, comparators = NULL, methods = NULL) {
  unique_values <- vapply(variables, function(col) {
    data.table::uniqueN(c(A[[col]], B[[col]]))
  }, integer(1))
  constant_vars <- names(unique_values)[unique_values == 1L]

  if (length(constant_vars) > 0L) {
    A[, (constant_vars) := NULL]
    B[, (constant_vars) := NULL]
    variables <- variables[!(variables %in% constant_vars)]

    if (!is.null(comparators)) {
      comparators[constant_vars] <- NULL
    }

    if (!is.null(methods)) {
      methods <- methods[!(names(methods) %in% constant_vars)]
    }
  }

  if (length(variables) == 0L) {
    stop("All key variables have only one unique value and cannot be used for record linkage.", call. = FALSE)
  }

  list(
    A = A,
    B = B,
    variables = variables,
    comparators = comparators,
    methods = methods,
    constant_vars = constant_vars
  )
}

#' @noRd
sanitize_true_matches <- function(true_matches, n_A, n_B, arg_name = "true_matches") {
  if (is.null(true_matches)) {
    return(NULL)
  }

  if (!(is.data.frame(true_matches) || is.data.table(true_matches))) {
    warning(sprintf("`%s` should be a data.frame or a data.table. Setting `%s` to `NULL`.", arg_name, arg_name))
    return(NULL)
  }

  true_matches <- data.table::copy(data.table::as.data.table(true_matches))

  tryCatch(
    {
      validate_match_pairs(true_matches, n_A, n_B, arg_name = arg_name)
      true_matches
    },
    error = function(e) {
      warning(sprintf("%s Setting `%s` to `NULL`.", conditionMessage(e), arg_name))
      NULL
    }
  )
}

#' @noRd
align_parameter_table <- function(params, variables) {
  params[match(variables, params[["variable"]])]
}

#' @noRd
fit_kliep_hurdle_model <- function(gamma_numerator, gamma_denominator, variable, controls_kliep) {
  gamma_num_dt <- positive_gamma_table(gamma_numerator, variable)
  gamma_den_dt <- positive_gamma_table(gamma_denominator, variable)

  if (NROW(gamma_num_dt) < 2L || NROW(gamma_den_dt) < 2L) {
    return(NULL)
  }

  do.call(
    densityratio::kliep,
    c(
      list(
        df_numerator = gamma_num_dt,
        df_denominator = gamma_den_dt
      ),
      controls_kliep
    )
  )
}

#' @noRd
fit_kliep_hurdle_models <- function(df_numerator, df_denominator, variables, controls_kliep) {
  ratio_kliep_list <- lapply(variables, function(variable) {
    fit_kliep_hurdle_model(
      gamma_numerator = df_numerator[[variable]],
      gamma_denominator = df_denominator[[variable]],
      variable = variable,
      controls_kliep = controls_kliep
    )
  })
  names(ratio_kliep_list) <- variables
  ratio_kliep_list
}

#' @noRd
missing_kliep_models <- function(ratio_kliep_list) {
  names(which(vapply(ratio_kliep_list, is.null, logical(1))))
}

#' @noRd
predict_kliep_positive <- function(model, gamma, variable) {
  kliep_pred <- rep(1, length(gamma))
  positive_idx <- which(gamma > 0)

  if (length(positive_idx) == 0L) {
    return(kliep_pred)
  }

  gamma_dt <- data.table::data.table(value = gamma[positive_idx])
  data.table::setnames(gamma_dt, "value", variable)
  kliep_pred[positive_idx] <- as.vector(stats::predict(model, gamma_dt))
  kliep_pred
}

#' @noRd
warn_kliep_issue <- function(context, variables, message) {
  warning(
    sprintf(
      "KLIEP issue in %s for variable%s %s: %s",
      context,
      if (length(variables) == 1L) "" else "s",
      paste(variables, collapse = ", "),
      message
    ),
    call. = FALSE
  )
}

#' @noRd
# Keep aligned parameter vectors to avoid repeated lookups inside the MEC loop.
bernoulli_product <- function(df, probs) {
  density_product <- rep(1, NROW(df))

  for (i in seq_along(probs)) {
    density_product <- density_product * stats::dbinom(
      x = df[[i]],
      size = 1,
      prob = probs[i]
    )
  }

  density_product
}

#' @noRd
bernoulli_ratio <- function(df, numerator_probs, denominator_probs) {
  bernoulli_product(df, numerator_probs) / bernoulli_product(df, denominator_probs)
}

#' @noRd
hurdle_gamma_product <- function(df, p_0, alpha, beta) {
  density_product <- rep(1, NROW(df))

  for (i in seq_along(p_0)) {
    density_product <- density_product * hurdle_gamma_density(
      x = df[[i]],
      p_0 = p_0[i],
      alpha = alpha[i],
      beta = beta[i]
    )
  }

  density_product
}

#' @noRd
hurdle_gamma_ratio <- function(df,
                               p_0_numerator,
                               alpha_numerator,
                               beta_numerator,
                               p_0_denominator,
                               alpha_denominator,
                               beta_denominator) {
  hurdle_gamma_product(df, p_0_numerator, alpha_numerator, beta_numerator) /
    hurdle_gamma_product(df, p_0_denominator, alpha_denominator, beta_denominator)
}

#' @noRd
kliep_hurdle_ratio <- function(df, variables, p_0_numerator, p_0_denominator, ratio_kliep_list) {
  ratio_list <- lapply(seq_along(variables), function(i) {
    variable <- variables[i]
    gamma_vec <- df[[variable]]

    if (!is.null(ratio_kliep_list[[variable]])) {
      kliep_pred <- predict_kliep_positive(ratio_kliep_list[[variable]], gamma_vec, variable)
      ifelse(gamma_vec == 0, p_0_numerator[i] / p_0_denominator[i], 1) *
        ifelse(gamma_vec > 0, (1 - p_0_numerator[i]) / (1 - p_0_denominator[i]) * kliep_pred, 1)
    } else {
      ifelse(gamma_vec == 0, p_0_numerator[i] / p_0_denominator[i], 1)
    }
  })

  Reduce(`*`, ratio_list)
}

#' @noRd
# Summarize the selected pairs and error rates without rebuilding large subsets in each step.
summarize_mec_selection <- function(a,
                                    b,
                                    ratio,
                                    g_est,
                                    n_M_est,
                                    duplicates_in_A = FALSE,
                                    set_construction = "size",
                                    target_rate = 0.03,
                                    tol = 0.005,
                                    max_iter = 50) {
  validate_choice(set_construction, c("size", "flr", "mmr"), "set_construction")

  order_idx <- order(ratio, decreasing = TRUE)
  ratio_sorted <- ratio[order_idx]
  g_est_sorted <- g_est[order_idx]
  g_est_cumsum <- cumsum(g_est_sorted)

  if (set_construction == "size") {
    selected_idx <- select_mec_indices(
      a = a,
      b = b,
      ratio = ratio,
      n_M = n_M_est,
      duplicates_in_A = duplicates_in_A
    )
    n_selected <- length(selected_idx)
    selected_g_est <- g_est[selected_idx]
    # Report the size-based rates for the returned classification set.
    flr_est <- if (n_selected == 0L) Inf else mean(1 - selected_g_est)
    mmr_est <- if (n_selected == 0L) 1 else max(0, min(1, 1 - sum(selected_g_est) / n_selected))

    return(list(
      selected_idx = selected_idx,
      flr_est = flr_est,
      mmr_est = mmr_est,
      iter = NULL
    ))
  }

  min_treshold <- min(ratio)
  max_treshold <- max(ratio)
  treshold <- (min_treshold + max_treshold) / 2
  iter <- 0

  while (iter < max_iter) {
    n_selected <- findInterval(-treshold, -ratio_sorted)

    if (set_construction == "flr") {
      current_rate <- if (n_selected == 0L) Inf else (n_selected - g_est_cumsum[n_selected]) / n_selected
    } else {
      current_rate <- if (n_selected == 0L || n_M_est <= 0) 1 else 1 - g_est_cumsum[n_selected] / n_M_est
    }

    if (abs(current_rate - target_rate) <= tol) {
      break
    } else if (current_rate < target_rate) {
      if (set_construction == "flr") {
        max_treshold <- treshold
      } else {
        min_treshold <- treshold
      }
    } else {
      if (set_construction == "flr") {
        min_treshold <- treshold
      } else {
        max_treshold <- treshold
      }
    }

    treshold <- (min_treshold + max_treshold) / 2
    iter <- iter + 1
  }

  n_selected <- findInterval(-treshold, -ratio_sorted)
  selected_idx <- if (n_selected == 0L) integer() else order_idx[seq_len(n_selected)]
  flr_est <- if (n_selected == 0L) Inf else (n_selected - g_est_cumsum[n_selected]) / n_selected
  mmr_est <- if (n_selected == 0L || n_M_est <= 0) 1 else 1 - g_est_cumsum[n_selected] / n_M_est

  list(
    selected_idx = selected_idx,
    flr_est = flr_est,
    mmr_est = mmr_est,
    iter = iter
  )
}

#' @noRd
# Rebuild the greedy MEC set from row indices instead of repeated table rbinding.
select_mec_indices <- function(a, b, ratio, n_M, duplicates_in_A = FALSE) {
  order_idx <- order(ratio, decreasing = TRUE)
  n_target <- round(n_M)

  if (n_target <= 0L || length(order_idx) == 0L) {
    return(integer())
  }

  selected_idx <- integer(length(order_idx))
  n_selected <- 0L
  used_a <- rep(FALSE, max(a))

  if (!duplicates_in_A) {
    used_b <- rep(FALSE, max(b))
  }

  for (idx in order_idx) {
    current_a <- a[idx]

    if (duplicates_in_A) {
      if (!used_a[current_a]) {
        n_selected <- n_selected + 1L
        selected_idx[n_selected] <- idx
        used_a[current_a] <- TRUE
      }
    } else {
      current_b <- b[idx]

      if (!used_a[current_a] && !used_b[current_b]) {
        n_selected <- n_selected + 1L
        selected_idx[n_selected] <- idx
        used_a[current_a] <- TRUE
        used_b[current_b] <- TRUE
      }
    }

    if (n_selected >= n_M) {
      break
    }
  }

  utils::head(selected_idx[seq_len(n_selected)], n_target)
}

#' @noRd
fixed_n_M <- function(n, ratio_gamma) {
  function(n_M) {
    sum(pmin(n_M * ratio_gamma / (n_M * (ratio_gamma - 1) + n), 1))
  }
}

#' @noRd
validate_blocking_controls <- function(controls_blocking) {
  if (!is.list(controls_blocking)) {
    stop("`controls_blocking` should be a list.")
  }

  protected_args <- intersect(names(controls_blocking), c("x", "y"))
  if (length(protected_args) > 0L) {
    stop("`controls_blocking` should not contain `x` or `y`.")
  }

  controls_blocking
}

#' @noRd
blocking_input_length <- function(x) {
  if (is.null(dim(x))) {
    length(x)
  } else {
    nrow(x)
  }
}

#' @noRd
build_blocking_input <- function(data, blocking_variables, blocking_sep) {
  if (!is.character(blocking_variables) || length(blocking_variables) == 0L) {
    stop("`blocking_variables` should be a non-empty character vector.")
  }

  if (!all(blocking_variables %in% names(data))) {
    stop("Not all `blocking_variables` are present in the input datasets.")
  }

  if (length(blocking_sep) != 1L || is.na(blocking_sep)) {
    stop("`blocking_sep` should be a single non-missing character value.")
  }

  if (length(blocking_variables) == 1L) {
    return(as.character(data[[blocking_variables]]))
  }

  blocking_data <- data[, blocking_variables, with = FALSE]
  do.call(paste, c(blocking_data, sep = blocking_sep))
}

#' @noRd
prepare_blocking_inputs <- function(A,
                                    B,
                                    blocking_x,
                                    blocking_y,
                                    blocking_variables,
                                    blocking_sep) {
  if (is.null(blocking_x) && is.null(blocking_y)) {
    blocking_x <- build_blocking_input(A, blocking_variables, blocking_sep)
    blocking_y <- build_blocking_input(B, blocking_variables, blocking_sep)
  } else if (is.null(blocking_x) || is.null(blocking_y)) {
    stop("`blocking_x` and `blocking_y` should both be supplied or both be `NULL`.")
  }

  if (blocking_input_length(blocking_x) != nrow(A)) {
    stop("`blocking_x` should have one row or element for each record in `A`.")
  }

  if (blocking_input_length(blocking_y) != nrow(B)) {
    stop("`blocking_y` should have one row or element for each record in `B`.")
  }

  list(x = blocking_x, y = blocking_y)
}

#' @noRd
run_blocking <- function(blocking_x, blocking_y, controls_blocking) {
  controls_blocking <- validate_blocking_controls(controls_blocking)

  do.call(
    blocking::blocking,
    c(
      list(x = blocking_x, y = blocking_y),
      controls_blocking
    )
  )
}

#' @noRd
extract_blocking_result_table <- function(blocking_result) {
  if (is.null(blocking_result[["result"]])) {
    stop("The object returned by `blocking::blocking()` should contain a `result` element.")
  }

  result <- data.table::as.data.table(blocking_result[["result"]])
  required_cols <- c("x", "y", "block")

  if (!all(required_cols %in% names(result))) {
    stop("`blocking::blocking()` result should contain columns x, y, and block.")
  }

  if (NROW(result) == 0L) {
    stop("`blocking::blocking()` did not return any candidate pairs.")
  }

  result[, required_cols, with = FALSE]
}

#' @noRd
reconstruct_block_summary <- function(blocking_table, n_A, n_B) {
  blocking_table <- data.table::copy(blocking_table)
  data.table::setnames(blocking_table, c("x", "y"), c("a", "b"))

  if (anyNA(blocking_table[["a"]]) || anyNA(blocking_table[["b"]]) ||
      anyNA(blocking_table[["block"]])) {
    stop("Blocking result cannot contain missing values in x, y, or block.")
  }

  if (!is.numeric(blocking_table[["a"]]) || !is.numeric(blocking_table[["b"]])) {
    stop("Blocking result columns x and y should contain numeric row indices.")
  }

  if (any(blocking_table[["a"]] != as.integer(blocking_table[["a"]])) ||
      any(blocking_table[["b"]] != as.integer(blocking_table[["b"]]))) {
    stop("Blocking result columns x and y should contain integer row indices.")
  }

  if (any(blocking_table[["a"]] < 1) || any(blocking_table[["a"]] > n_A) ||
      any(blocking_table[["b"]] < 1) || any(blocking_table[["b"]] > n_B)) {
    stop("Blocking result contains row indices outside the input datasets.")
  }

  x_blocks <- unique(blocking_table[, .(a, block)])
  y_blocks <- unique(blocking_table[, .(b, block)])

  if (anyDuplicated(x_blocks[["a"]]) > 0L ||
      anyDuplicated(y_blocks[["b"]]) > 0L) {
    stop("Final blocks should be disjoint on both input files.")
  }

  A_by_block <- x_blocks[order(block, a), .(A = list(a), n_A = .N), by = block]
  B_by_block <- y_blocks[order(block, b), .(B = list(b), n_B = .N), by = block]
  block_summary <- merge(A_by_block, B_by_block, by = "block", all = FALSE)
  data.table::setorder(block_summary, block)
  block_summary[, pair_count := n_A * n_B]
  block_summary[, nonmatches_min := pair_count - pmin(n_A, n_B)]

  excluded_records <- list(
    A = setdiff(seq_len(n_A), x_blocks[["a"]]),
    B = setdiff(seq_len(n_B), y_blocks[["b"]])
  )

  list(
    block_summary = block_summary,
    excluded_records = excluded_records
  )
}

#' @noRd
make_block_pair_table <- function(block_summary, blocks = NULL) {
  if (is.null(blocks)) {
    selected_blocks <- block_summary[["block"]]
  } else {
    selected_blocks <- blocks
  }

  block_summary <- block_summary[block %in% selected_blocks]
  pair_list <- lapply(seq_len(NROW(block_summary)), function(i) {
    pairs <- data.table::CJ(
      a = block_summary[["A"]][[i]],
      b = block_summary[["B"]][[i]]
    )
    data.table::setkey(pairs, NULL)
    pairs[, block := block_summary[["block"]][i]]
    pairs[]
  })

  data.table::rbindlist(pair_list, use.names = TRUE)
}

#' @noRd
sample_with_optional_seed <- function(n, seed = NULL) {
  if (is.null(seed)) {
    return(sample.int(n))
  }

  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
  }

  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })

  set.seed(seed)
  sample.int(n)
}

#' @noRd
validate_nonmatch_sample_size <- function(nonmatch_sample_size, n_A, n_B) {
  n_full <- as.numeric(n_A) * as.numeric(n_B)

  if (is.null(nonmatch_sample_size)) {
    return(n_full)
  }

  if (length(nonmatch_sample_size) != 1L || is.na(nonmatch_sample_size) ||
      nonmatch_sample_size <= 0 || nonmatch_sample_size != floor(nonmatch_sample_size)) {
    stop("`nonmatch_sample_size` should be `NULL` or a positive integer.")
  }

  if (nonmatch_sample_size > n_full) {
    stop("`nonmatch_sample_size` cannot be larger than `nrow(A) * nrow(B)`.")
  }

  as.numeric(nonmatch_sample_size)
}

#' @noRd
sample_cartesian_pair_ids <- function(n_pairs, sample_size, seed = NULL) {
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
  }

  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })

  if (!is.null(seed)) {
    set.seed(seed)
  }

  sample.int(n_pairs, sample_size, replace = FALSE)
}

#' @noRd
cartesian_pair_ids_to_pairs <- function(pair_ids, n_B) {
  data.table::data.table(
    a = ((pair_ids - 1) %/% n_B) + 1,
    b = ((pair_ids - 1) %% n_B) + 1
  )
}

#' @noRd
update_nonmatch_summary <- function(summary, Omega, b_vars, cpar_vars) {
  if (length(b_vars) > 0L) {
    summary$b_sum <- summary$b_sum + colSums(as.matrix(Omega[, b_vars, with = FALSE]))
  }

  if (length(cpar_vars) > 0L) {
    for (var in cpar_vars) {
      gamma <- Omega[[var]]
      positive_gamma <- gamma[gamma > 0]
      summary$cpar_zero[var] <- summary$cpar_zero[var] + sum(gamma == 0)
      summary$cpar_positive_n[var] <- summary$cpar_positive_n[var] + length(positive_gamma)
      summary$cpar_positive_sum[var] <- summary$cpar_positive_sum[var] + sum(positive_gamma)
      summary$cpar_log_positive_sum[var] <- summary$cpar_log_positive_sum[var] + sum(log(positive_gamma))
    }
  }

  summary
}

#' @noRd
estimate_nonmatch_parameters <- function(A,
                                         B,
                                         variables,
                                         comparators,
                                         methods,
                                         nonmatch_sample_size,
                                         nonmatch_sampling_seed = NULL,
                                         controls_nleqslv = list(),
                                         chunk_size = 50000L) {
  n_A <- nrow(A)
  n_B <- nrow(B)
  n_full <- as.numeric(n_A) * as.numeric(n_B)
  actual_sample_size <- validate_nonmatch_sample_size(nonmatch_sample_size, n_A, n_B)
  method_variables <- extract_method_variables(methods, include_hit_miss = FALSE)
  b_vars <- method_variables$b_vars
  cpar_vars <- method_variables$cpar_vars

  summary <- list(
    b_sum = stats::setNames(numeric(length(b_vars)), b_vars),
    cpar_zero = stats::setNames(numeric(length(cpar_vars)), cpar_vars),
    cpar_positive_n = stats::setNames(numeric(length(cpar_vars)), cpar_vars),
    cpar_positive_sum = stats::setNames(numeric(length(cpar_vars)), cpar_vars),
    cpar_log_positive_sum = stats::setNames(numeric(length(cpar_vars)), cpar_vars)
  )

  sampled_ids <- NULL
  if (actual_sample_size < n_full) {
    sampled_ids <- sample_cartesian_pair_ids(
      n_pairs = n_full,
      sample_size = actual_sample_size,
      seed = nonmatch_sampling_seed
    )
  }

  start <- 1
  while (start <= actual_sample_size) {
    end <- min(start + chunk_size - 1, actual_sample_size)
    pair_ids <- if (is.null(sampled_ids)) {
      start:end
    } else {
      sampled_ids[start:end]
    }
    pairs <- cartesian_pair_ids_to_pairs(pair_ids, n_B)
    vectors <- comparison_vectors(
      A = A,
      B = B,
      variables = variables,
      comparators = comparators,
      pairs = pairs
    )
    summary <- update_nonmatch_summary(
      summary = summary,
      Omega = vectors$Omega,
      b_vars = b_vars,
      cpar_vars = cpar_vars
    )
    start <- end + 1
  }

  b_params <- NULL
  if (length(b_vars) > 0L) {
    b_params <- data.table::data.table(
      variable = b_vars,
      eta = summary$b_sum / actual_sample_size
    )
  }

  cpar_params <- NULL
  if (length(cpar_vars) > 0L) {
    gamma_plus_U <- summary$cpar_positive_sum / summary$cpar_positive_n
    modified_nleqslv <- purrr::partial(nleqslv::nleqslv, control = controls_nleqslv)
    alpha_U <- alpha_formula_summary(
      n_positive = summary$cpar_positive_n,
      positive_sum = summary$cpar_positive_sum,
      log_positive_sum = summary$cpar_log_positive_sum,
      fun = modified_nleqslv
    )
    beta_U <- alpha_U / gamma_plus_U

    cpar_params <- data.table::data.table(
      variable = cpar_vars,
      p_0_U = summary$cpar_zero / actual_sample_size,
      alpha_U = alpha_U,
      beta_U = beta_U
    )
  }

  list(
    sample_size = actual_sample_size,
    binary = b_params,
    continuous_parametric = cpar_params
  )
}

#' @noRd
select_training_blocks <- function(block_summary,
                                   min_training_pairs,
                                   min_training_nonmatches,
                                   block_sampling_seed = NULL) {
  if (is.null(min_training_pairs) && is.null(min_training_nonmatches)) {
    training_blocks <- data.table::copy(block_summary)
    training_blocks[, cumulative_pairs := cumsum(pair_count)]
    training_blocks[, cumulative_nonmatches_min := cumsum(nonmatches_min)]
    return(list(
      training_rule = "all_blocks",
      training_blocks = training_blocks
    ))
  }

  if (is.null(min_training_pairs) || is.null(min_training_nonmatches)) {
    stop("`min_training_pairs` and `min_training_nonmatches` should both be supplied or both be `NULL`.")
  }

  if (length(min_training_pairs) != 1L || is.na(min_training_pairs) ||
      min_training_pairs <= 0) {
    stop("`min_training_pairs` should be a positive number.")
  }

  if (length(min_training_nonmatches) != 1L || is.na(min_training_nonmatches) ||
      min_training_nonmatches <= 0) {
    stop("`min_training_nonmatches` should be a positive number.")
  }

  permutation <- sample_with_optional_seed(NROW(block_summary), block_sampling_seed)
  training_blocks <- data.table::copy(block_summary[permutation])
  training_blocks[, cumulative_pairs := cumsum(pair_count)]
  training_blocks[, cumulative_nonmatches_min := cumsum(nonmatches_min)]
  hit_idx <- which(
    training_blocks[["cumulative_pairs"]] >= min_training_pairs &
      training_blocks[["cumulative_nonmatches_min"]] >= min_training_nonmatches
  )

  if (length(hit_idx) == 0L) {
    stop("Training sample is infeasible under the current blocking configuration and thresholds.")
  }

  training_blocks <- training_blocks[seq_len(hit_idx[1L])]
  list(
    training_rule = "threshold_sampling",
    training_blocks = training_blocks
  )
}

#' @noRd
exact_match_pairs <- function(A, B, variables) {
  merge(
    A[, c("a", variables), with = FALSE],
    B[, c("b", variables), with = FALSE],
    by = variables
  )[, c("a", "b"), with = FALSE]
}

#' @noRd
score_mec_ratio <- function(Omega, model) {
  data.table::set(Omega, j = "ratio", value = 1)

  if (!is.null(model$b_vars)) {
    b_vars <- model$b_vars
    b_params <- align_parameter_table(model$b_params, b_vars)
    Omega_b <- Omega[, b_vars, with = FALSE]
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * bernoulli_ratio(Omega_b, b_params$theta, b_params$eta)
    )
  }

  if (!is.null(model$cpar_vars)) {
    cpar_vars <- model$cpar_vars
    cpar_params <- data.table::copy(align_parameter_table(model$cpar_params, cpar_vars))
    if ("p_0_U" %in% names(cpar_params)) {
      p_0_denominator <- cpar_params$p_0_U
      alpha_denominator <- cpar_params$alpha_U
      beta_denominator <- cpar_params$beta_U
    } else {
      p_0_denominator <- cpar_params$p_0_Omega
      alpha_denominator <- cpar_params$alpha_Omega
      beta_denominator <- cpar_params$beta_Omega
    }
    Omega_cpar <- Omega[, cpar_vars, with = FALSE]
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * hurdle_gamma_ratio(
        Omega_cpar,
        cpar_params$p_0_M,
        cpar_params$alpha_M,
        cpar_params$beta_M,
        p_0_denominator,
        alpha_denominator,
        beta_denominator
      )
    )
  }

  if (!is.null(model$cnonpar_vars)) {
    cnonpar_vars <- model$cnonpar_vars
    Omega_cnonpar <- Omega[, cnonpar_vars, with = FALSE]

    if (!is.null(model$ratio_kliep_list)) {
      ratio_kliep_list <- model$ratio_kliep_list
      p_0_M_cnonpar <- model$cnonpar_params[["p_0_M_cnonpar"]]
      p_0_U_cnonpar <- model$cnonpar_params[["p_0_U_cnonpar"]]
      names(p_0_M_cnonpar) <- cnonpar_vars
      names(p_0_U_cnonpar) <- cnonpar_vars

      ratio_kliep <- kliep_hurdle_ratio(
        Omega_cnonpar,
        cnonpar_vars,
        p_0_M_cnonpar,
        p_0_U_cnonpar,
        ratio_kliep_list
      )
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * ratio_kliep)
    } else if (!is.null(model$ratio_kliep)) {
      data.table::set(
        Omega,
        j = "ratio",
        value = Omega[["ratio"]] * stats::predict(model$ratio_kliep, Omega_cnonpar)
      )
    }
  }

  if (!is.null(model$hm_vars)) {
    hm_vars <- model$hm_vars
    hm_params <- align_parameter_table(model$hm_params, hm_vars)
    Omega_hm <- Omega[, hm_vars, with = FALSE]
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * bernoulli_ratio(Omega_hm, hm_params$theta, hm_params$eta)
    )
  }

  Omega
}

#' @noRd
estimate_local_mec_size <- function(ratio,
                                    n_pairs,
                                    n_A,
                                    n_B,
                                    fixed_method,
                                    prob_ratio,
                                    prob_est = NULL) {
  n_M_max <- min(n_A, n_B)

  if (prob_ratio == "1") {
    if (is.null(prob_est)) {
      stop("`prob_est` is required when `prob_ratio == \"1\"`.")
    }
    g_est <- pmin(prob_est * ratio, 1)
    g_est[!is.finite(g_est)] <- 0
    n_M_est <- min(sum(g_est), n_M_max)
    n_M_est <- max(n_M_est, 0)
    n_M_est <- round(n_M_est)
    g_est <- pmin(prob_est * ratio, 1)
    g_est[!is.finite(g_est)] <- 0
  } else {
    fun_n_M <- fixed_n_M(n = n_pairs, ratio_gamma = ratio)
    n_M_original <- tryCatch(
      FixedPoint::FixedPoint(
        Function = fun_n_M,
        Inputs = n_M_max,
        Method = fixed_method
      )$FixedPoint,
      error = function(e) NA_real_
    )
    if (length(n_M_original) != 1L || !is.finite(n_M_original)) {
      n_M_original <- n_M_max
    }
    n_M_est <- min(n_M_original, n_M_max)
    n_M_est <- max(n_M_est, 0)
    n_M_est <- round(n_M_est)
    g_est <- pmin(n_M_est * ratio / (n_M_est * (ratio - 1) + n_pairs), 1)
    g_est[!is.finite(g_est)] <- 0
  }

  list(n_M_est = n_M_est, g_est = g_est)
}

#' @noRd
# Singleton blocks have no one-to-one conflict, so classify the only pair by
# its fitted density ratio.
select_singleton_mec_index <- function(ratio) {
  ratio <- ratio[1L]

  if (!is.na(ratio) && ratio >= 1) {
    return(1L)
  }

  integer()
}

#' @noRd
blocking_diagnostics <- function(true_matches, block_pairs, n_full_pairs) {
  true_matches <- data.table::copy(true_matches)
  data.table::setkey(true_matches, a, b)
  block_pairs_key <- data.table::copy(block_pairs[, c("a", "b"), with = FALSE])
  data.table::setkey(block_pairs_key, a, b)

  preserved <- true_matches[block_pairs_key, nomatch = 0L]
  n_true <- NROW(true_matches)
  n_preserved <- NROW(preserved)

  data.table::data.table(
    true_matches = n_true,
    preserved_matches = n_preserved,
    lost_matches = n_true - n_preserved,
    blocking_recall = if (n_true == 0L) NA_real_ else n_preserved / n_true,
    blocking_fnr = if (n_true == 0L) NA_real_ else 1 - n_preserved / n_true,
    blocked_pairs = NROW(block_pairs),
    full_pairs = n_full_pairs
  )
}

#' @noRd
mec_selection_diagnostics <- function(true_matches, block_pairs, M_est) {
  true_matches <- data.table::copy(data.table::as.data.table(true_matches))
  block_pairs <- data.table::copy(data.table::as.data.table(block_pairs))[, .(a, b)]
  M_pairs <- data.table::copy(data.table::as.data.table(M_est))[, .(a, b)]

  candidate_matches <- NROW(block_pairs[true_matches, on = .(a, b), nomatch = 0L])
  selected_candidate_matches <- NROW(M_pairs[true_matches, on = .(a, b), nomatch = 0L])
  selected_pairs <- NROW(M_pairs)
  true_match_count <- NROW(true_matches)
  blocking_lost_matches <- true_match_count - candidate_matches
  mec_missed_candidate_matches <- candidate_matches - selected_candidate_matches
  false_links <- selected_pairs - selected_candidate_matches

  data.table::data.table(
    true_matches = true_match_count,
    candidate_matches = candidate_matches,
    blocking_lost_matches = blocking_lost_matches,
    selected_candidate_matches = selected_candidate_matches,
    mec_missed_candidate_matches = mec_missed_candidate_matches,
    mec_candidate_recall = if (candidate_matches == 0L) NA_real_ else selected_candidate_matches / candidate_matches,
    mec_candidate_fnr = if (candidate_matches == 0L) NA_real_ else mec_missed_candidate_matches / candidate_matches,
    selected_pairs = selected_pairs,
    false_links = false_links,
    empirical_flr = if (selected_pairs == 0L) NA_real_ else false_links / selected_pairs
  )
}
