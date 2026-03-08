#' @import data.table
#' @importFrom stats rbinom
#' @importFrom stats rgamma
#' @importFrom stats dbinom
#' @importFrom utils head
#'
#' @title Parametric Bootstrap for Standard Error Estimation in MEC
#'
#' @author Adam Struzik
#'
#' @description
#' **WORK IN PROGRESS**
#'
#' *This function is currently under development.*
#'
#' @param mec_result An object of class `mec_rec_lin`.
#' @param B A number of bootstrap iterations.
#'
#' @export
est_se_bootstrap <- function(mec_result,
                             B = 100) {

  # TODO add `mec_rec_lin` class validation

  # TODO add information about no support for cnonpar and hm

  n <- mec_result$n
  n_M_original <- mec_result$n_M_est
  n_M_est <- round(n_M_original)

  b_vars <- mec_result$b_vars
  cpar_vars <- mec_result$cpar_vars

  b_params <- mec_result$b_params
  cpar_params <- mec_result$cpar_params

  delta <- mec_result$delta
  eps <- mec_result$eps
  controls_nleqslv <- mec_result$controls_nleqslv

  n_M_boot <- numeric(B)

  for (b in seq_len(B)) {

    g <- sample(c(rep(1, n_M_est), rep(0, n - n_M_est)))
    Omega <- data.table(g = g)

    if (length(b_vars) > 0) {
      for (var in b_vars) {

        theta <- as.numeric(b_params[b_params[["variable"]] == var, "theta"])
        eta <- as.numeric(b_params[b_params[["variable"]] == var, "eta"])

        prob_vec <- fifelse(g == 1, theta, eta)
        set(Omega, j = var, value = rbinom(n, 1, prob_vec))

      }
    }

    if (length(cpar_vars) > 0) {
      for(var in cpar_vars) {

        p_0_M <- as.numeric(cpar_params[cpar_params[["variable"]] == var, "p_0_M"])
        alpha_M <- as.numeric(cpar_params[cpar_params[["variable"]] == var, "alpha_M"])
        beta_M <- as.numeric(cpar_params[cpar_params[["variable"]] == var, "beta_M"])
        p_0_U <- as.numeric(cpar_params[cpar_params[["variable"]] == var, "p_0_U"])
        alpha_U <- as.numeric(cpar_params[cpar_params[["variable"]] == var, "alpha_U"])
        beta_U <- as.numeric(cpar_params[cpar_params[["variable"]] == var, "beta_U"])

        p_0_vec <- fifelse(g == 1, p_0_M, p_0_U)
        alpha_vec <- fifelse(g == 1, alpha_M, alpha_U)
        beta_vec <- fifelse(g == 1, beta_M, beta_U)

        is_positive <- rbinom(n, 1, 1 - p_0_vec)
        gamma_values <- rgamma(n, shape = alpha_vec, rate = beta_vec)
        final_values <- fifelse(is_positive == 1, gamma_values, 0)
        set(Omega, j = var, value = final_values)

      }
    }

    n_M_boot[b] <- mec_relaxed(
      Omega = Omega,
      n = n,
      n_M_est = n_M_est,
      b_vars = b_vars,
      cpar_vars = cpar_vars,
      b_params_init = b_params,
      cpar_params_init = cpar_params,
      controls_nleqslv = controls_nleqslv,
      delta = delta,
      eps = eps
    )

  }

  se_n_M <- sqrt(1 / (B - 1) * sum((n_M_boot - n_M_original)^2))

  return(se_n_M)

}

#' @noRd
mec_relaxed <- function(Omega,
                        n,
                        n_M_est,
                        b_vars,
                        cpar_vars,
                        b_params_init,
                        cpar_params_init,
                        controls_nleqslv,
                        delta,
                        eps) {

  b_params <- copy(b_params_init)
  cpar_params <- copy(cpar_params_init)

  set(Omega, j = "ratio", value = 1)

  if (length(b_vars) > 0) {

    Omega_b <- Omega[, b_vars, with = FALSE]

    b_params$eta <- binary_formula(Omega_b)

    b_numerator_list <- lapply(b_vars, function(col) {
      dbinom(x = Omega_b[[col]],
             size = 1,
             prob = as.numeric(b_params[b_params[["variable"]] == col, "theta"]))
    })
    b_numerator <- Reduce(`*`, b_numerator_list)
    b_denominator_list <- lapply(b_vars, function(col) {
      dbinom(x = Omega_b[[col]],
             size = 1,
             prob = as.numeric(b_params[b_params[["variable"]] == col, "eta"]))
    })
    b_denominator <- Reduce(`*`, b_denominator_list)

    set(Omega, j = "ratio", value = Omega[["ratio"]] * b_numerator / b_denominator)
    set(Omega, j = "b_denominator", value = b_denominator)

  }

  if (length(cpar_vars) > 0) {

    Omega_cpar <- Omega[, cpar_vars, with = FALSE]

    cpar_params$p_0_U <- p_0_formula(Omega_cpar)
    gamma_plus_U <- gamma_plus_formula(Omega_cpar)
    modified_nleqslv <- purrr::partial(nleqslv::nleqslv, control = controls_nleqslv)
    alpha_U <- alpha_formula(Omega_cpar, modified_nleqslv)
    cpar_params$alpha_U <- alpha_U
    cpar_params$beta_U <- alpha_U / gamma_plus_U

    cpar_numerator_list <- lapply(cpar_vars, function(col) {
      hurdle_gamma_density(x = Omega_cpar[[col]],
                           p_0 = as.numeric(cpar_params[cpar_params[["variable"]] == col, "p_0_M"]),
                           alpha = as.numeric(cpar_params[cpar_params[["variable"]] == col, "alpha_M"]),
                           beta = as.numeric(cpar_params[cpar_params[["variable"]] == col, "beta_M"]))
    })
    cpar_numerator <- Reduce(`*`, cpar_numerator_list)
    cpar_denominator_list <- lapply(cpar_vars, function(col) {
      hurdle_gamma_density(x = Omega_cpar[[col]],
                           p_0 = as.numeric(cpar_params[cpar_params[["variable"]] == col, "p_0_U"]),
                           alpha = as.numeric(cpar_params[cpar_params[["variable"]] == col, "alpha_U"]),
                           beta = as.numeric(cpar_params[cpar_params[["variable"]] == col, "beta_U"]))
    })
    cpar_denominator <- Reduce(`*`, cpar_denominator_list)

    set(Omega, j = "ratio", value = Omega[["ratio"]] * cpar_numerator / cpar_denominator)
    set(Omega, j = "cpar_denominator", value = cpar_denominator)

  }

  iter <- 1

  repeat {

    if (iter == 1) {
      g_est <- pmin(n_M_est * Omega[["ratio"]] / (n_M_est * (Omega[["ratio"]] - 1) + n), 1)
      n_M_old <- n_M_est
    } else {
      g_est <- pmin(NROW(M) * Omega[["ratio"]] / (NROW(M) * (Omega[["ratio"]] - 1) + n), 1)
      n_M_old <- n_M
    }

    n_M <- sum(g_est)
    set(Omega, j = "g_est", value = g_est)

    Omega <- Omega[order(-get("ratio")), ]
    M <- head(Omega, round(n_M))

    if (NROW(M) == 0) {
      break
    }

    if (iter >= 2) {

      old_params <- c()
      params <- c()
      if (length(b_vars) > 0) {
        old_params <- c(old_params, theta_b_old)
        params <- c(params, theta_b)
      }
      if (length(cpar_vars) > 0) {
        old_params <- c(old_params, p_0_M_old, alpha_M_old, beta_M_old)
        params <- c(params, p_0_M, alpha_M, beta_M)
      }
      if ((abs(n_M_old - n_M) < delta) || norm(old_params - params, type = "2") < eps) {
        break
      }

    }

    set(Omega, j = "ratio", value = 1)

    if (length(b_vars) > 0) {

      Omega_b <- Omega[, b_vars, with = FALSE]
      M_b <- M[, b_vars, with = FALSE]

      theta_b_old <- b_params$theta
      theta_b <- binary_formula(M_b)

      b_params$theta <- theta_b

      b_numerator_list <- lapply(b_vars, function(col) {
        stats::dbinom(x = Omega[[col]],
                      size = 1,
                      prob = as.numeric(b_params[b_params[["variable"]] == col, "theta"]))
      })
      b_numerator <- Reduce(`*`, b_numerator_list)
      set(Omega, j = "ratio", value = Omega[["ratio"]] * b_numerator / Omega[["b_denominator"]])

    }

    if (length(cpar_vars) > 0) {

      Omega_cpar <- Omega[, cpar_vars, with = FALSE]
      M_cpar <- M[, cpar_vars, with = FALSE]

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

      cpar_numerator_list <- lapply(cpar_vars, function(col) {
        hurdle_gamma_density(x = Omega[[col]],
                             p_0 = as.numeric(cpar_params[cpar_params[["variable"]] == col, "p_0_M"]),
                             alpha = as.numeric(cpar_params[cpar_params[["variable"]] == col, "alpha_M"]),
                             beta = as.numeric(cpar_params[cpar_params[["variable"]] == col, "beta_M"]))
      })
      cpar_numerator <- Reduce(`*`, cpar_numerator_list)

      set(Omega, j = "ratio", value = Omega[["ratio"]] * cpar_numerator / Omega[["cpar_denominator"]])

    }

    iter <- iter + 1

  }

  return(n_M)

}
