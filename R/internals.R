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
      methods[constant_vars] <- NULL
    }
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
        ifelse(gamma_vec > 0, (1 - p_0_numerator[i]) * (1 - p_0_denominator[i]) * kliep_pred, 1)
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
    selected_g_est <- g_est[selected_idx]
    flr_est <- if (length(selected_idx) == 0L) Inf else mean(1 - selected_g_est)
    mmr_est <- if (n_M_est <= 0) 1 else 1 - sum(selected_g_est) / n_M_est

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
