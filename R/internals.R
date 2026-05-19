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

#' @noRd
fit_mec_unsupervised_omega <- function(A,
                                       B,
                                       variables,
                                       comparators,
                                       methods,
                                       Omega,
                                       exact_matches,
                                       start_params = NULL,
                                       nonmatch_params = NULL,
                                       nonpar_hurdle = TRUE,
                                       delta = 0.5,
                                       eps = 0.05,
                                       max_iter_em = 10,
                                       tol_em = 1,
                                       controls_nleqslv = list(),
                                       controls_kliep = control_kliep(),
                                       context = "mec_blocking()") {
  if (!is.null(start_params)) {
    stopifnot("`start_params` should be a list." =
                is.list(start_params))
  }
  if (!is.null(nonmatch_params)) {
    stopifnot("`nonmatch_params` should be a list." =
                is.list(nonmatch_params))
  }

  method_variables <- extract_method_variables(methods, include_hit_miss = TRUE)
  b_vars <- method_variables$b_vars
  cpar_vars <- method_variables$cpar_vars
  cnonpar_vars <- method_variables$cnonpar_vars
  hm_vars <- method_variables$hm_vars

  exact_match_idx <- data.table(
    omega_idx = seq_len(NROW(Omega)),
    a = Omega[["a"]],
    b = Omega[["b"]]
  )[exact_matches[, c("a", "b"), with = FALSE], on = c("a", "b"), nomatch = 0L][["omega_idx"]]

  if (length(exact_match_idx) == 0L) {
    stop("The MEC training space contains no exact agreements on the key variables.")
  }

  n <- NROW(Omega)
  exact_match_mask <- rep(FALSE, n)
  exact_match_mask[exact_match_idx] <- TRUE
  M_idx <- exact_match_idx
  n_M <- length(M_idx)
  all_idx <- seq_len(n)
  selection_mask <- rep(FALSE, n)

  b_params <- NULL
  cpar_params <- NULL
  cnonpar_params <- NULL
  hm_params <- NULL
  ratio_kliep <- NULL
  ratio_kliep_list <- NULL
  kliep_warning_emitted <- FALSE

  warn_kliep_once <- function(variables, message) {
    if (!kliep_warning_emitted) {
      warn_kliep_issue(context, variables, message)
      kliep_warning_emitted <<- TRUE
    }
  }

  data.table::set(Omega, j = "ratio", value = 1)

  if (is.null(start_params)) {
    start_params <- list()

    if (length(b_vars) > 0L) {
      start_params[["binary"]] <- data.table(
        variable = b_vars,
        theta = runif(length(b_vars), min = 0.9)
      )
    }

    if (length(cpar_vars) > 0L) {
      start_params[["continuous_parametric"]] <- data.table(
        variable = cpar_vars,
        p_0_M = runif(length(cpar_vars), min = 0.8, max = 0.9),
        alpha_M = runif(length(cpar_vars), min = 0.1, max = 1),
        beta_M = runif(length(cpar_vars), min = 10, max = 20)
      )
    }

    if (length(hm_vars) > 0L) {
      start_params[["hit_miss"]] <- data.table(
        variable = hm_vars,
        theta = runif(length(hm_vars), min = 0.9)
      )
    }
  }

  if (length(b_vars) > 0L) {
    b_params <- align_parameter_table(start_params$binary, b_vars)
    Omega_b <- Omega[, b_vars, with = FALSE]
    if (!is.null(nonmatch_params) && !is.null(nonmatch_params$binary)) {
      eta_b <- align_parameter_table(nonmatch_params$binary, b_vars)$eta
    } else {
      eta_b <- binary_formula(Omega_b)
    }
    b_params$eta <- eta_b
    b_denominator <- bernoulli_product(Omega_b, b_params$eta)
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * bernoulli_ratio(Omega_b, b_params$theta, b_params$eta)
    )
    data.table::set(Omega, j = "b_denominator", value = b_denominator)
  }

  if (length(cpar_vars) > 0L) {
    cpar_params <- align_parameter_table(start_params$continuous_parametric, cpar_vars)
    Omega_cpar <- Omega[, cpar_vars, with = FALSE]
    modified_nleqslv <- purrr::partial(nleqslv::nleqslv, control = controls_nleqslv)
    if (!is.null(nonmatch_params) && !is.null(nonmatch_params$continuous_parametric)) {
      cpar_nonmatch_params <- align_parameter_table(nonmatch_params$continuous_parametric, cpar_vars)
      p_0_U <- cpar_nonmatch_params$p_0_U
      alpha_U <- cpar_nonmatch_params$alpha_U
      beta_U <- cpar_nonmatch_params$beta_U
    } else {
      p_0_U <- p_0_formula(Omega_cpar)
      gamma_plus_U <- gamma_plus_formula(Omega_cpar)
      alpha_U <- alpha_formula(Omega_cpar, modified_nleqslv)
      beta_U <- alpha_U / gamma_plus_U
    }
    cpar_params$p_0_U <- p_0_U
    cpar_params$alpha_U <- alpha_U
    cpar_params$beta_U <- beta_U
    cpar_denominator <- hurdle_gamma_product(
      Omega_cpar,
      cpar_params$p_0_U,
      cpar_params$alpha_U,
      cpar_params$beta_U
    )
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * hurdle_gamma_ratio(
        Omega_cpar,
        cpar_params$p_0_M,
        cpar_params$alpha_M,
        cpar_params$beta_M,
        cpar_params$p_0_U,
        cpar_params$alpha_U,
        cpar_params$beta_U
      )
    )
    data.table::set(Omega, j = "cpar_denominator", value = cpar_denominator)
  }

  if (length(cnonpar_vars) > 0L) {
    Omega_cnonpar <- Omega[, cnonpar_vars, with = FALSE]
    n_exact_matches <- sum(exact_match_mask)

    if (nonpar_hurdle) {
      p_0_M_cnonpar <- runif(length(cnonpar_vars), min = 0.5)
      p_0_U_cnonpar <- p_0_formula(Omega_cnonpar)
      names(p_0_M_cnonpar) <- cnonpar_vars
      cnonpar_params <- data.table(
        variable = cnonpar_vars,
        p_0_M_cnonpar = p_0_M_cnonpar,
        p_0_U_cnonpar = p_0_U_cnonpar
      )

      ratio_temp <- lapply(cnonpar_vars, function(x) {
        r <- numeric(n)
        r[exact_match_mask] <- stats::runif(n_exact_matches, min = 5, max = 10)
        r[!exact_match_mask] <- stats::runif(n - n_exact_matches, min = 0.1, max = 1)
        r
      })
      names(ratio_temp) <- cnonpar_vars
      cnonpar_ratio_list <- lapply(cnonpar_vars, function(x) {
        gamma_vec <- Omega_cnonpar[[x]]
        ifelse(gamma_vec == 0, p_0_M_cnonpar[x] / p_0_U_cnonpar[x], 1) *
          ifelse(gamma_vec > 0, (1 - p_0_M_cnonpar[x]) / (1 - p_0_U_cnonpar[x]) * ratio_temp[[x]], 1)
      })
      ratio_kliep <- Reduce(`*`, cnonpar_ratio_list)
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * ratio_kliep)
    } else {
      ratio_temp <- as.numeric(Omega$ratio)
      ratio_temp[exact_match_mask] <- ratio_temp[exact_match_mask] * stats::runif(
        n_exact_matches,
        min = 5,
        max = 10
      )
      ratio_temp[!exact_match_mask] <- ratio_temp[!exact_match_mask] * stats::runif(
        n - n_exact_matches,
        min = 0.1,
        max = 5
      )
      data.table::set(Omega, j = "ratio", value = ratio_temp)
    }
  }

  if (length(hm_vars) > 0L) {
    hm_params <- align_parameter_table(start_params$hit_miss, hm_vars)
    Omega_hm <- Omega[, hm_vars, with = FALSE]
    eta_hm <- binary_formula(Omega_hm)
    hm_params$eta <- eta_hm
    hm_denominator <- bernoulli_product(Omega_hm, hm_params$eta)
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * bernoulli_ratio(Omega_hm, hm_params$theta, hm_params$eta)
    )
    data.table::set(Omega, j = "hm_denominator", value = hm_denominator)

    values_list_m <- lapply(hm_vars, function(x) {
      name <- substr(x, 7, nchar(x))
      values <- unique(c(A[[name]], B[[name]]))
      m_est <- sapply(values, function(y) sum(A[[name]] == y) / NROW(A))
      data.table::data.table(
        "value" = values,
        "m_est" = m_est
      )
    })
    names(values_list_m) <- hm_vars
  }

  iter <- 1

  repeat {
    g_est <- pmin(length(M_idx) * Omega$ratio / (length(M_idx) * (Omega$ratio - 1) + n), 1)
    n_M_old <- n_M
    n_M <- sum(g_est)

    if (n_M > min(NROW(A), NROW(B))) {
      n_M <- min(NROW(A), NROW(B))
    }

    data.table::set(Omega, j = "g_est", value = g_est)
    M_idx <- select_mec_indices(
      a = Omega[["a"]],
      b = Omega[["b"]],
      ratio = Omega[["ratio"]],
      n_M = n_M,
      duplicates_in_A = FALSE
    )

    if (length(M_idx) == 0L) {
      break
    }

    selection_mask[M_idx] <- TRUE
    U_idx <- all_idx[!selection_mask]
    selection_mask[M_idx] <- FALSE

    if (iter >= 2) {
      old_params <- c()
      params <- c()
      if (length(b_vars) > 0L) {
        old_params <- c(old_params, theta_b_old)
        params <- c(params, theta_b)
      }
      if (length(cpar_vars) > 0L) {
        old_params <- c(old_params, p_0_M_old, alpha_M_old, beta_M_old)
        params <- c(params, p_0_M, alpha_M, beta_M)
      }
      if (length(hm_vars) > 0L) {
        old_params <- c(old_params, theta_hm_old)
        params <- c(params, theta_hm)
      }

      if (length(cnonpar_vars) == 0L) {
        if ((abs(n_M_old - n_M) < delta) || norm(old_params - params, type = "2") < eps) {
          break
        }
      } else if ((abs(n_M_old - n_M) < delta)) {
        break
      }
    }

    data.table::set(Omega, j = "ratio", value = 1)

    if (length(b_vars) > 0L) {
      M_b <- Omega_b[M_idx]
      theta_b_old <- b_params$theta
      theta_b <- binary_formula(M_b)
      b_params$theta <- theta_b
      b_numerator <- bernoulli_product(Omega_b, b_params$theta)
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * b_numerator / Omega[["b_denominator"]])
    }

    if (length(cpar_vars) > 0L) {
      M_cpar <- Omega_cpar[M_idx]
      p_0_M_old <- cpar_params$p_0_M
      alpha_M_old <- cpar_params$alpha_M
      beta_M_old <- cpar_params$beta_M
      p_0_M <- p_0_formula(M_cpar)
      gamma_plus_M <- gamma_plus_formula(M_cpar)
      alpha_M <- alpha_formula(M_cpar, modified_nleqslv)
      beta_M <- alpha_M / gamma_plus_M
      beta_M[is.nan(beta_M)] <- beta_M_old[is.nan(beta_M)]
      cpar_params$p_0_M <- p_0_M
      cpar_params$alpha_M <- alpha_M
      cpar_params$beta_M <- beta_M
      cpar_numerator <- hurdle_gamma_product(
        Omega_cpar,
        cpar_params$p_0_M,
        cpar_params$alpha_M,
        cpar_params$beta_M
      )
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * cpar_numerator / Omega[["cpar_denominator"]])
    }

    if (length(cnonpar_vars) > 0L) {
      M_cnonpar <- Omega_cnonpar[M_idx]
      U_cnonpar <- Omega_cnonpar[U_idx]

      if (nonpar_hurdle) {
        p_0_M_cnonpar <- p_0_formula(M_cnonpar)
        p_0_U_cnonpar <- p_0_formula(U_cnonpar)
        cnonpar_params <- data.table(
          variable = cnonpar_vars,
          p_0_M_cnonpar = p_0_M_cnonpar,
          p_0_U_cnonpar = p_0_U_cnonpar
        )
        ratio_kliep_old <- ratio_kliep

        tryCatch({
          ratio_kliep_models <- fit_kliep_hurdle_models(
            df_numerator = M_cnonpar,
            df_denominator = U_cnonpar,
            variables = cnonpar_vars,
            controls_kliep = controls_kliep
          )
          missing_models <- missing_kliep_models(ratio_kliep_models)

          if (length(missing_models) < length(cnonpar_vars)) {
            if (length(missing_models) > 0L) {
              warn_kliep_once(missing_models, "using only the hurdle mass term for those variables in the current iteration.")
            }

            ratio_kliep_list <- ratio_kliep_models
            ratio_kliep <- kliep_hurdle_ratio(
              df = Omega_cnonpar,
              variables = cnonpar_vars,
              p_0_numerator = p_0_M_cnonpar,
              p_0_denominator = p_0_U_cnonpar,
              ratio_kliep_list = ratio_kliep_models
            )
          } else {
            ratio_kliep <- ratio_kliep_old
            warn_kliep_once(cnonpar_vars, "could not be fitted in the current iteration; using the previous ratio estimate.")
          }
          data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * ratio_kliep)
        },
          error = function(e) {
            warn_kliep_once(cnonpar_vars, paste("failed with error:", e$message))
          })
      } else {
        ratio_kliep <- do.call(
          densityratio::kliep,
          c(list(
            df_numerator = M_cnonpar,
            df_denominator = U_cnonpar
          ),
          controls_kliep)
        )
        data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * stats::predict(ratio_kliep, Omega_cnonpar))
      }
    }

    if (length(hm_vars) > 0L) {
      M_hm <- Omega_hm[M_idx]
      theta_hm_old <- hm_params$theta
      theta_hm <- binary_formula(M_hm)
      hm_params$theta <- theta_hm

      p_est <- n_M / max(NROW(A), NROW(B))
      u_est_list <- lapply(hm_vars, function(x) {
        u_est <- runif(NROW(values_list_m[[x]]), 0, 1)
        u_est <- u_est / sum(u_est)
      })
      names(u_est_list) <- hm_vars

      m_bk_prod <- sapply(seq_len(NROW(B)), function(x) {
        m_bk_list <- lapply(hm_vars, function(y) {
          var_name <- substr(y, 7, nchar(y))
          d_value <- (B[[var_name]])[x]
          as.vector((values_list_m[[y]])[(values_list_m[[y]])[["value"]] == d_value, ][["m_est"]])
        })
        names(m_bk_list) <- hm_vars
        Reduce(`*`, m_bk_list)
      })

      iter_em <- 1
      while (iter_em <= max_iter_em + 1) {
        if (iter_em >= 2) {
          delta_b_old <- em_params[["delta_b"]]
        }

        delta_m_u <- lapply(seq_len(NROW(B)), function(x) {
          u_bk_list <- lapply(hm_vars, function(y) {
            var_name <- substr(y, 7, nchar(y))
            values <- cbind(values_list_m[[y]], u_est_list[[y]])
            colnames(values) <- c("value", "m_est", "u_est")
            d_value <- (B[[var_name]])[x]
            as.vector(values[values[["value"]] == d_value, ][["u_est"]])
          })
          names(u_bk_list) <- hm_vars
          u_bk_prod <- Reduce(`*`, u_bk_list)
          data.table::data.table(
            delta_b = p_est * m_bk_prod[x] / (p_est * m_bk_prod[x] + (1 - p_est) * u_bk_prod),
            m_bk_prod = m_bk_prod[x],
            u_bk_prod = u_bk_prod
          )
        })

        em_params <- data.table::rbindlist(delta_m_u)

        if (iter_em >= 3) {
          log_lik_old <- log_lik
        }

        if (iter_em >= 2) {
          log_lik <- sum(
            ifelse(delta_b_old == 0,
                   0,
                   delta_b_old * log(p_est * em_params[["m_bk_prod"]]))
          ) +
            sum(ifelse(em_params[["u_bk_prod"]] == 0,
                       0,
                       (1 - delta_b_old) * log((1 - p_est) * em_params[["u_bk_prod"]]))
                )
        }

        if (iter_em >= 3) {
          if (abs(log_lik - log_lik_old) <= tol_em) break
        }

        u_est_list <- lapply(hm_vars, function(x) {
          var_name <- substr(x, 7, nchar(x))
          z_bk <- B[[var_name]]

          sapply((values_list_m[[x]])[["value"]], function(y) {
            sum((1 - em_params[["delta_b"]]) * (z_bk == y)) / sum(1 - em_params[["delta_b"]])
          })
        })
        names(u_est_list) <- hm_vars

        iter_em <- iter_em + 1
      }

      hm_numerator <- bernoulli_product(Omega_hm, hm_params$theta)
      eta_hm <- sapply(hm_vars, function(col) {
        ((1 - p_est) * sum(u_est_list[[col]] * (values_list_m[[col]])[["m_est"]]) +
           p_est * (1 - 1 / NROW(A)) * sum(((values_list_m[[col]])[["m_est"]]) ^ 2)) /
          (1 - p_est / NROW(A))
      })
      hm_params$eta <- eta_hm
      hm_denominator <- bernoulli_product(Omega_hm, hm_params$eta)
      data.table::set(Omega, j = "hm_denominator", value = hm_denominator)
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * hm_numerator / Omega[["hm_denominator"]])
    }

    iter <- iter + 1
  }

  n_M_est <- n_M
  selection_summary <- summarize_mec_selection(
    a = Omega[["a"]],
    b = Omega[["b"]],
    ratio = Omega[["ratio"]],
    g_est = Omega[["g_est"]],
    n_M_est = n_M_est,
    duplicates_in_A = FALSE,
    set_construction = "size"
  )
  M <- Omega[selection_summary$selected_idx]

  model <- list(
    b_vars = if (length(b_vars) == 0L) NULL else b_vars,
    cpar_vars = if (length(cpar_vars) == 0L) NULL else cpar_vars,
    cnonpar_vars = if (length(cnonpar_vars) == 0L) NULL else cnonpar_vars,
    hm_vars = if (length(hm_vars) == 0L) NULL else hm_vars,
    b_params = if (is.null(b_params)) NULL else b_params,
    cpar_params = if (is.null(cpar_params)) NULL else cpar_params,
    cnonpar_params = if (is.null(cnonpar_params)) NULL else cnonpar_params,
    hm_params = if (is.null(hm_params)) NULL else hm_params,
    ratio_kliep = if (is.null(ratio_kliep) || is.numeric(ratio_kliep)) NULL else ratio_kliep,
    ratio_kliep_list = if (is.null(ratio_kliep_list)) NULL else ratio_kliep_list,
    variables = variables,
    comparators = comparators,
    methods = methods,
    n_M_est = n_M_est,
    prob_est = n_M_est / n
  )

  list(
    model = model,
    Omega = Omega,
    M_est = M[, c("a", "b", "ratio"), with = FALSE],
    n_M_est = n_M_est,
    flr_est = selection_summary$flr_est,
    mmr_est = selection_summary$mmr_est
  )
}

#' @noRd
validate_mec_blocking_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha < 0 || alpha >= 1) {
    stop("`alpha` should be a single numeric value in [0, 1).", call. = FALSE)
  }

  invisible(alpha)
}

#' @noRd
initial_inverted_match_count <- function(nu) {
  as.integer(max(0L, nu))
}

#' @noRd
prepare_inverted_start_params <- function(start_params, cpar_vars) {
  if (!is.null(start_params)) {
    stopifnot("`start_params` should be a list." =
                is.list(start_params))
    return(start_params)
  }

  start_params <- list()
  if (length(cpar_vars) > 0L) {
    start_params[["continuous_parametric"]] <- data.table(
      variable = cpar_vars,
      p_0_M = stats::runif(length(cpar_vars), min = 0.8, max = 0.9),
      alpha_M = stats::runif(length(cpar_vars), min = 0.1, max = 1),
      beta_M = stats::runif(length(cpar_vars), min = 10, max = 20)
    )
  }

  start_params
}

#' @noRd
get_hurdle_gamma_fallback <- function(fallback_params, variables, side) {
  if (is.null(fallback_params)) {
    return(NULL)
  }

  alpha_col <- paste0("alpha_", side)
  beta_col <- paste0("beta_", side)
  if (!all(c("variable", alpha_col, beta_col) %in% names(fallback_params))) {
    return(NULL)
  }

  fallback_params <- align_parameter_table(fallback_params, variables)
  list(
    alpha = fallback_params[[alpha_col]],
    beta = fallback_params[[beta_col]]
  )
}

#' @noRd
estimate_hurdle_gamma_params <- function(df,
                                         variables,
                                         side,
                                         fallback_params = NULL,
                                         controls_nleqslv = list(),
                                         context = "mec_blocking()") {
  if (length(variables) == 0L) {
    return(NULL)
  }

  p_0 <- p_0_formula(df)
  alpha <- stats::setNames(rep(NA_real_, length(variables)), variables)
  beta <- stats::setNames(rep(NA_real_, length(variables)), variables)
  fallback <- get_hurdle_gamma_fallback(fallback_params, variables, side)
  modified_nleqslv <- purrr::partial(nleqslv::nleqslv, control = controls_nleqslv)

  for (i in seq_along(variables)) {
    variable <- variables[i]
    gamma <- df[[variable]]
    positive_gamma <- gamma[gamma > 0]

    if (length(positive_gamma) >= 2L) {
      estimated_alpha <- tryCatch(
        alpha_formula(df[, variable, with = FALSE], modified_nleqslv)[1L],
        error = function(e) NA_real_
      )
      estimated_beta <- estimated_alpha / mean(positive_gamma)

      if (is.finite(estimated_alpha) && is.finite(estimated_beta) &&
          estimated_alpha > 0 && estimated_beta > 0) {
        alpha[i] <- estimated_alpha
        beta[i] <- estimated_beta
        next
      }
    }

    has_fallback <- !is.null(fallback) &&
      is.finite(fallback$alpha[i]) && is.finite(fallback$beta[i]) &&
      fallback$alpha[i] > 0 && fallback$beta[i] > 0

    if (!has_fallback) {
      stop(
        sprintf(
          "%s cannot estimate %s-side Gamma parameters for `%s`; at least two positive continuous comparisons or finite fallback parameters are required.",
          context,
          if (side == "M") "match" else "nonmatch",
          variable
        ),
        call. = FALSE
      )
    }

    alpha[i] <- fallback$alpha[i]
    beta[i] <- fallback$beta[i]
  }

  params <- data.table(variable = variables)
  data.table::set(params, j = paste0("p_0_", side), value = as.numeric(p_0))
  data.table::set(params, j = paste0("alpha_", side), value = as.numeric(alpha))
  data.table::set(params, j = paste0("beta_", side), value = as.numeric(beta))
  params
}

#' @noRd
estimate_inverted_match_parameters <- function(Omega,
                                               M_idx,
                                               b_vars,
                                               cpar_vars,
                                               previous_params,
                                               start_params,
                                               controls_nleqslv,
                                               context = "mec_blocking()") {
  if (length(M_idx) == 0L) {
    M_idx <- seq_len(NROW(Omega))
  }

  match_Omega <- Omega[M_idx]
  cpar_fallback <- start_params$continuous_parametric
  if (!is.null(previous_params$continuous_parametric)) {
    cpar_fallback <- previous_params$continuous_parametric
  }

  b_params <- NULL
  if (length(b_vars) > 0L) {
    b_params <- data.table(
      variable = b_vars,
      theta = binary_formula(match_Omega[, b_vars, with = FALSE])
    )
  }

  cpar_params <- NULL
  if (length(cpar_vars) > 0L) {
    cpar_params <- estimate_hurdle_gamma_params(
      df = match_Omega[, cpar_vars, with = FALSE],
      variables = cpar_vars,
      side = "M",
      fallback_params = cpar_fallback,
      controls_nleqslv = controls_nleqslv,
      context = context
    )
  }

  list(binary = b_params, continuous_parametric = cpar_params)
}

#' @noRd
estimate_inverted_nonmatch_parameters <- function(Omega,
                                                  U_idx,
                                                  b_vars,
                                                  cpar_vars,
                                                  previous_params,
                                                  controls_nleqslv,
                                                  context = "mec_blocking()") {
  if (length(U_idx) == 0L) {
    stop(sprintf("%s cannot estimate nonmatch parameters from an empty candidate-pair complement.", context),
         call. = FALSE)
  }

  b_params <- NULL
  if (length(b_vars) > 0L) {
    b_params <- data.table(
      variable = b_vars,
      eta = binary_formula(Omega[U_idx, b_vars, with = FALSE])
    )
  }

  cpar_params <- NULL
  if (length(cpar_vars) > 0L) {
    cpar_params <- estimate_hurdle_gamma_params(
      df = Omega[U_idx, cpar_vars, with = FALSE],
      variables = cpar_vars,
      side = "U",
      fallback_params = previous_params$continuous_parametric,
      controls_nleqslv = controls_nleqslv,
      context = context
    )
  }

  list(binary = b_params, continuous_parametric = cpar_params)
}

#' @noRd
combine_inverted_parameters <- function(match_params, nonmatch_params, b_vars, cpar_vars) {
  b_params <- NULL
  if (length(b_vars) > 0L) {
    b_params <- data.table(
      variable = b_vars,
      theta = match_params$binary$theta,
      eta = nonmatch_params$binary$eta
    )
  }

  cpar_params <- NULL
  if (length(cpar_vars) > 0L) {
    cpar_params <- data.table::copy(match_params$continuous_parametric)
    data.table::set(cpar_params, j = "p_0_U", value = nonmatch_params$continuous_parametric$p_0_U)
    data.table::set(cpar_params, j = "alpha_U", value = nonmatch_params$continuous_parametric$alpha_U)
    data.table::set(cpar_params, j = "beta_U", value = nonmatch_params$continuous_parametric$beta_U)
  }

  list(binary = b_params, continuous_parametric = cpar_params)
}

#' @noRd
inverted_nonmatch_param_vector <- function(nonmatch_params) {
  params <- numeric()

  if (!is.null(nonmatch_params$binary)) {
    params <- c(params, nonmatch_params$binary$eta)
  }

  if (!is.null(nonmatch_params$continuous_parametric)) {
    params <- c(
      params,
      nonmatch_params$continuous_parametric$p_0_U,
      nonmatch_params$continuous_parametric$alpha_U,
      nonmatch_params$continuous_parametric$beta_U
    )
  }

  params
}

#' @noRd
inverted_match_param_vector <- function(match_params) {
  params <- numeric()

  if (!is.null(match_params$binary)) {
    params <- c(params, match_params$binary$theta)
  }

  if (!is.null(match_params$continuous_parametric)) {
    params <- c(
      params,
      match_params$continuous_parametric$p_0_M,
      match_params$continuous_parametric$alpha_M,
      match_params$continuous_parametric$beta_M
    )
  }

  params
}

#' @noRd
inverted_parameter_vector <- function(match_params, nonmatch_params) {
  c(
    inverted_match_param_vector(match_params),
    inverted_nonmatch_param_vector(nonmatch_params)
  )
}

#' @noRd
score_inverted_mec_ratio <- function(Omega, b_vars, cpar_vars, b_params, cpar_params) {
  data.table::set(Omega, j = "ratio", value = 1)

  if (length(b_vars) > 0L) {
    Omega_b <- Omega[, b_vars, with = FALSE]
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * bernoulli_ratio(Omega_b, b_params$eta, b_params$theta)
    )
  }

  if (length(cpar_vars) > 0L) {
    Omega_cpar <- Omega[, cpar_vars, with = FALSE]
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * hurdle_gamma_ratio(
        Omega_cpar,
        cpar_params$p_0_U,
        cpar_params$alpha_U,
        cpar_params$beta_U,
        cpar_params$p_0_M,
        cpar_params$alpha_M,
        cpar_params$beta_M
      )
    )
  }

  ratio <- Omega[["ratio"]]
  ratio[is.na(ratio) | ratio < 0] <- Inf
  data.table::set(Omega, j = "ratio", value = ratio)
  Omega
}

#' @noRd
blocking_disagreement_norm <- function(Omega, b_vars, cpar_vars) {
  components <- list()

  if (length(b_vars) > 0L) {
    components[["binary"]] <- 1 - as.matrix(Omega[, b_vars, with = FALSE])
  }

  if (length(cpar_vars) > 0L) {
    components[["continuous_parametric"]] <- as.matrix(Omega[, cpar_vars, with = FALSE])
  }

  disagreement <- do.call(cbind, components)
  sqrt(rowSums(disagreement ^ 2))
}

#' @noRd
select_inverted_mec_indices <- function(a, b, block, ratio, n_M) {
  n_target <- round(n_M)

  if (n_target <= 0L || length(ratio) == 0L) {
    return(integer())
  }

  ratio_order <- ratio
  ratio_order[is.na(ratio_order) | ratio_order < 0] <- Inf
  order_idx <- order(ratio_order, a, b, block)
  selected_idx <- integer(length(order_idx))
  n_selected <- 0L
  used_a <- rep(FALSE, max(a))
  used_b <- rep(FALSE, max(b))

  for (idx in order_idx) {
    current_a <- a[idx]
    current_b <- b[idx]

    if (!used_a[current_a] && !used_b[current_b]) {
      n_selected <- n_selected + 1L
      selected_idx[n_selected] <- idx
      used_a[current_a] <- TRUE
      used_b[current_b] <- TRUE
    }

    if (n_selected >= n_target) {
      break
    }
  }

  if (n_selected == 0L) {
    return(integer())
  }

  selected_idx[seq_len(n_selected)]
}

#' @noRd
estimate_inverted_q <- function(ratio, n_U_current, N) {
  denominator <- n_U_current * (ratio - 1) + N
  q_est <- n_U_current * ratio / denominator
  q_est[is.infinite(ratio) & ratio > 0] <- 1
  q_est[ratio == 0 & denominator > 0] <- 0
  q_est <- pmin(q_est, 1)
  q_est <- pmax(q_est, 0)
  q_est[!is.finite(q_est)] <- 1
  q_est
}

#' @noRd
has_valid_u_cpar_fallback <- function(previous_params, cpar_vars) {
  if (length(cpar_vars) == 0L) {
    return(logical())
  }

  if (is.null(previous_params)) {
    return(stats::setNames(rep(FALSE, length(cpar_vars)), cpar_vars))
  }

  cpar_params <- previous_params$continuous_parametric
  if (is.null(cpar_params) ||
      !all(c("variable", "alpha_U", "beta_U") %in% names(cpar_params))) {
    return(stats::setNames(rep(FALSE, length(cpar_vars)), cpar_vars))
  }

  cpar_params <- align_parameter_table(cpar_params, cpar_vars)
  valid <- is.finite(cpar_params[["alpha_U"]]) &
    is.finite(cpar_params[["beta_U"]]) &
    cpar_params[["alpha_U"]] > 0 &
    cpar_params[["beta_U"]] > 0
  stats::setNames(valid, cpar_vars)
}

#' @noRd
has_sufficient_u_fit_sample <- function(Omega,
                                        U_fit_idx,
                                        cpar_vars,
                                        previous_params) {
  if (length(U_fit_idx) == 0L) {
    return(FALSE)
  }

  if (length(cpar_vars) == 0L) {
    return(TRUE)
  }

  valid_fallback <- has_valid_u_cpar_fallback(previous_params, cpar_vars)
  for (variable in cpar_vars) {
    gamma <- Omega[[variable]][U_fit_idx]
    n_positive <- sum(is.finite(gamma) & gamma > 0)
    if (n_positive < 2L && !isTRUE(valid_fallback[[variable]])) {
      return(FALSE)
    }
  }

  TRUE
}

#' @noRd
make_u_fit_diagnostic <- function(iter,
                                  n_U_current,
                                  n_U_base,
                                  alpha,
                                  requested_n_drop,
                                  requested_n_keep,
                                  actual_n_drop,
                                  n_U_fit,
                                  alpha_applied,
                                  reason) {
  data.table(
    iter = as.integer(iter),
    n_U_current = as.integer(n_U_current),
    n_U_base = as.integer(n_U_base),
    alpha = as.numeric(alpha),
    requested_n_drop = as.integer(requested_n_drop),
    requested_n_keep = as.integer(requested_n_keep),
    actual_n_drop = as.integer(actual_n_drop),
    actual_drop_fraction = if (n_U_base == 0L) NA_real_ else actual_n_drop / n_U_base,
    n_U_fit = as.integer(n_U_fit),
    alpha_applied = as.logical(alpha_applied),
    reason = reason
  )
}

#' @noRd
rank_u_fit_candidates <- function(Omega, U_fit_idx) {
  score_name <- if ("q_est" %in% names(Omega)) "q_est" else "ratio"
  score <- Omega[[score_name]][U_fit_idx]
  if (all(is.na(score)) && "ratio" %in% names(Omega)) {
    score <- Omega[["ratio"]][U_fit_idx]
  }
  score[is.na(score)] <- -Inf
  U_fit_idx[order(-score, Omega[["a"]][U_fit_idx], Omega[["b"]][U_fit_idx], Omega[["block"]][U_fit_idx])]
}

#' @noRd
select_inverted_u_fit_indices <- function(Omega,
                                          U_idx,
                                          U_base_idx,
                                          iter,
                                          alpha,
                                          cpar_vars,
                                          previous_params) {
  n_U_current <- length(U_idx)
  n_U_base <- length(U_base_idx)
  requested_n_drop <- floor(alpha * n_U_current)
  requested_n_keep <- n_U_current - requested_n_drop

  if (iter == 1L) {
    return(list(
      U_fit_idx = U_base_idx,
      diagnostic = make_u_fit_diagnostic(
        iter = iter,
        n_U_current = n_U_current,
        n_U_base = n_U_base,
        alpha = alpha,
        requested_n_drop = 0L,
        requested_n_keep = n_U_current,
        actual_n_drop = 0L,
        n_U_fit = n_U_base,
        alpha_applied = FALSE,
        reason = "first_u_fit_full"
      )
    ))
  }

  if (alpha == 0 || requested_n_drop == 0L) {
    reason <- if (alpha == 0) "alpha_zero" else "requested_drop_zero"
    return(list(
      U_fit_idx = U_base_idx,
      diagnostic = make_u_fit_diagnostic(
        iter = iter,
        n_U_current = n_U_current,
        n_U_base = n_U_base,
        alpha = alpha,
        requested_n_drop = requested_n_drop,
        requested_n_keep = requested_n_keep,
        actual_n_drop = 0L,
        n_U_fit = n_U_base,
        alpha_applied = FALSE,
        reason = reason
      )
    ))
  }

  if (n_U_base <= requested_n_keep) {
    return(list(
      U_fit_idx = U_base_idx,
      diagnostic = make_u_fit_diagnostic(
        iter = iter,
        n_U_current = n_U_current,
        n_U_base = n_U_base,
        alpha = alpha,
        requested_n_drop = requested_n_drop,
        requested_n_keep = requested_n_keep,
        actual_n_drop = 0L,
        n_U_fit = n_U_base,
        alpha_applied = FALSE,
        reason = "base_smaller_than_requested_keep"
      )
    ))
  }

  ranked_idx <- rank_u_fit_candidates(Omega, U_base_idx)
  retained_idx <- ranked_idx[seq_len(requested_n_keep)]
  if (!has_sufficient_u_fit_sample(
    Omega = Omega,
    U_fit_idx = retained_idx,
    cpar_vars = cpar_vars,
    previous_params = previous_params
  )) {
    return(list(
      U_fit_idx = U_base_idx,
      diagnostic = make_u_fit_diagnostic(
        iter = iter,
        n_U_current = n_U_current,
        n_U_base = n_U_base,
        alpha = alpha,
        requested_n_drop = requested_n_drop,
        requested_n_keep = requested_n_keep,
        actual_n_drop = 0L,
        n_U_fit = n_U_base,
        alpha_applied = FALSE,
        reason = "minimum_sample_full_base"
      )
    ))
  }

  actual_n_drop <- n_U_base - length(retained_idx)
  list(
    U_fit_idx = retained_idx,
    diagnostic = make_u_fit_diagnostic(
      iter = iter,
      n_U_current = n_U_current,
      n_U_base = n_U_base,
      alpha = alpha,
      requested_n_drop = requested_n_drop,
      requested_n_keep = requested_n_keep,
      actual_n_drop = actual_n_drop,
      n_U_fit = length(retained_idx),
      alpha_applied = actual_n_drop > 0L,
      reason = "alpha_reliability_drop"
    )
  )
}

#' @noRd
fit_mec_blocking_inverted_omega <- function(A,
                                            B,
                                            variables,
                                            comparators,
                                            methods,
                                            Omega,
                                            start_params = NULL,
                                            delta = 0.5,
                                            eps = 0.05,
                                            controls_nleqslv = list(),
                                            n_U_min,
                                            nu,
                                            alpha = 0,
                                            context = "mec_blocking()") {
  method_variables <- extract_method_variables(methods, include_hit_miss = FALSE)
  b_vars <- method_variables$b_vars
  cpar_vars <- method_variables$cpar_vars
  N <- NROW(Omega)
  all_idx <- seq_len(N)
  start_params <- prepare_inverted_start_params(start_params, cpar_vars)
  max_iter_inverted <- 1000L

  init_disagreement <- blocking_disagreement_norm(Omega, b_vars, cpar_vars)
  data.table::set(Omega, j = "init_disagreement", value = init_disagreement)
  n_M_init_target <- initial_inverted_match_count(nu)
  M_idx <- select_inverted_mec_indices(
    a = Omega[["a"]],
    b = Omega[["b"]],
    block = Omega[["block"]],
    ratio = Omega[["init_disagreement"]],
    n_M = n_M_init_target
  )
  U_idx <- setdiff(all_idx, M_idx)
  n_M_init <- length(M_idx)
  n_U_init <- length(U_idx)

  if (length(U_idx) == 0L) {
    if (N != nu) {
      stop(sprintf("%s initialized an empty nonmatch complement before reaching the structural one-to-one bound.",
                   context),
           call. = FALSE)
    }

    match_params <- estimate_inverted_match_parameters(
      Omega = Omega,
      M_idx = M_idx,
      b_vars = b_vars,
      cpar_vars = cpar_vars,
      previous_params = list(),
      start_params = start_params,
      controls_nleqslv = controls_nleqslv,
      context = context
    )
    data.table::set(Omega, j = "ratio", value = 0)
    data.table::set(Omega, j = "q_est", value = 0)
    M_est <- Omega[M_idx, c("a", "b", "block", "ratio"), with = FALSE]
    u_fit_diagnostics <- make_u_fit_diagnostic(
      iter = 0L,
      n_U_current = 0L,
      n_U_base = 0L,
      alpha = alpha,
      requested_n_drop = 0L,
      requested_n_keep = 0L,
      actual_n_drop = 0L,
      n_U_fit = 0L,
      alpha_applied = FALSE,
      reason = "structural_no_nonmatch_complement"
    )
    model <- list(
      b_vars = if (length(b_vars) == 0L) NULL else b_vars,
      cpar_vars = if (length(cpar_vars) == 0L) NULL else cpar_vars,
      b_params = match_params$binary,
      cpar_params = match_params$continuous_parametric,
      nonmatch_params = NULL,
      variables = variables,
      comparators = comparators,
      methods = methods,
      n_M_est = NROW(M_est),
      n_U_est = 0L,
      n_U_min = n_U_min,
      nu = nu,
      n_M_init = n_M_init,
      n_U_init = n_U_init,
      candidate_pair_count = N,
      prob_est = if (N == 0L) NA_real_ else NROW(M_est) / N,
      ratio_orientation = "u_over_m",
      alpha = alpha,
      n_U_fit = 0L,
      u_fit_diagnostics = u_fit_diagnostics,
      iter = 0L,
      convergence_reason = "structural_no_nonmatch_complement"
    )

    return(list(
      model = model,
      Omega = Omega,
      M_est = M_est,
      n_M_est = NROW(M_est),
      n_U_est = 0L,
      n_U_min = n_U_min,
      nu = nu,
      n_M_init = n_M_init,
      n_U_init = n_U_init,
      candidate_pair_count = N,
      alpha = alpha,
      n_U_fit = 0L,
      u_fit_diagnostics = u_fit_diagnostics,
      iter = 0L,
      convergence_reason = "structural_no_nonmatch_complement"
    ))
  }

  iter <- 1L
  n_U_old <- length(U_idx)
  previous_match_params <- NULL
  previous_nonmatch_params <- NULL
  previous_param_vector <- NULL
  convergence_reason <- "max_iter"
  n_U_fit <- length(U_idx)
  u_fit_diagnostics <- data.table()

  repeat {
    match_params <- estimate_inverted_match_parameters(
      Omega = Omega,
      M_idx = M_idx,
      b_vars = b_vars,
      cpar_vars = cpar_vars,
      previous_params = if (is.null(previous_match_params)) list() else previous_match_params,
      start_params = start_params,
      controls_nleqslv = controls_nleqslv,
      context = context
    )
    u_fit_selection <- select_inverted_u_fit_indices(
      Omega = Omega,
      U_idx = U_idx,
      U_base_idx = U_idx,
      iter = iter,
      alpha = alpha,
      cpar_vars = cpar_vars,
      previous_params = previous_nonmatch_params
    )
    U_fit_idx <- u_fit_selection$U_fit_idx
    u_fit_diagnostics <- rbindlist(
      list(u_fit_diagnostics, u_fit_selection$diagnostic),
      use.names = TRUE
    )
    n_U_fit <- length(U_fit_idx)
    nonmatch_params <- estimate_inverted_nonmatch_parameters(
      Omega = Omega,
      U_idx = U_fit_idx,
      b_vars = b_vars,
      cpar_vars = cpar_vars,
      previous_params = previous_nonmatch_params,
      controls_nleqslv = controls_nleqslv,
      context = context
    )
    combined_params <- combine_inverted_parameters(
      match_params = match_params,
      nonmatch_params = nonmatch_params,
      b_vars = b_vars,
      cpar_vars = cpar_vars
    )
    Omega <- score_inverted_mec_ratio(
      Omega = Omega,
      b_vars = b_vars,
      cpar_vars = cpar_vars,
      b_params = combined_params$binary,
      cpar_params = combined_params$continuous_parametric
    )

    q_est <- estimate_inverted_q(
      ratio = Omega[["ratio"]],
      n_U_current = length(U_idx),
      N = N
    )
    data.table::set(Omega, j = "q_est", value = q_est)
    n_U_raw <- sum(q_est)
    n_U_est <- max(n_U_min, min(N, round(n_U_raw)))
    n_M_est <- N - n_U_est
    M_idx_new <- select_inverted_mec_indices(
      a = Omega[["a"]],
      b = Omega[["b"]],
      block = Omega[["block"]],
      ratio = Omega[["ratio"]],
      n_M = n_M_est
    )
    U_idx_new <- setdiff(all_idx, M_idx_new)
    param_vector <- inverted_parameter_vector(match_params, nonmatch_params)
    can_check_convergence <- iter >= 2L

    if (can_check_convergence && abs(n_U_est - n_U_old) < delta) {
      convergence_reason <- "n_U_delta"
    } else if (can_check_convergence && identical(sort(M_idx_new), sort(M_idx))) {
      convergence_reason <- "match_set_unchanged"
    } else if (can_check_convergence &&
               !is.null(previous_param_vector) &&
               norm(previous_param_vector - param_vector, type = "2") < eps) {
      convergence_reason <- "nonmatch_parameter_eps"
    } else if (iter >= max_iter_inverted) {
      warning("mec_blocking() reached the maximum number of inverted MEC iterations.")
      convergence_reason <- "max_iter"
    } else if (length(U_idx_new) == 0L) {
      convergence_reason <- "structural_no_nonmatch_complement"
    } else if (length(M_idx_new) == 0L) {
      convergence_reason <- "empty_match_set"
    } else {
      previous_match_params <- match_params
      previous_nonmatch_params <- nonmatch_params
      previous_param_vector <- param_vector
      n_U_old <- n_U_est
      M_idx <- M_idx_new
      U_idx <- U_idx_new
      iter <- iter + 1L
      next
    }

    M_idx <- M_idx_new
    U_idx <- U_idx_new
    previous_match_params <- match_params
    previous_nonmatch_params <- nonmatch_params
    previous_param_vector <- param_vector
    break
  }

  M_est <- Omega[M_idx, c("a", "b", "block", "ratio"), with = FALSE]
  n_M_selected <- NROW(M_est)
  n_U_selected <- N - n_M_selected
  model <- list(
    b_vars = if (length(b_vars) == 0L) NULL else b_vars,
    cpar_vars = if (length(cpar_vars) == 0L) NULL else cpar_vars,
    b_params = combined_params$binary,
    cpar_params = combined_params$continuous_parametric,
    nonmatch_params = previous_nonmatch_params,
    variables = variables,
    comparators = comparators,
    methods = methods,
    n_M_est = n_M_selected,
    n_U_est = n_U_selected,
    n_U_min = n_U_min,
    nu = nu,
    n_M_init = n_M_init,
    n_U_init = n_U_init,
    candidate_pair_count = N,
    prob_est = n_M_selected / N,
    ratio_orientation = "u_over_m",
    alpha = alpha,
    n_U_fit = n_U_fit,
    u_fit_diagnostics = u_fit_diagnostics,
    iter = iter,
    convergence_reason = convergence_reason
  )

  list(
    model = model,
    Omega = Omega,
    M_est = M_est,
    n_M_est = n_M_selected,
    n_U_est = n_U_selected,
    n_U_min = n_U_min,
    nu = nu,
    n_M_init = n_M_init,
    n_U_init = n_U_init,
    candidate_pair_count = N,
    alpha = alpha,
    n_U_fit = n_U_fit,
    u_fit_diagnostics = u_fit_diagnostics,
    iter = iter,
    convergence_reason = convergence_reason
  )
}
