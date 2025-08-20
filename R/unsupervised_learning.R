#' @import data.table
#' @importFrom stats runif
#' @importFrom stats predict
#' @importFrom purrr partial
#' @importFrom utils head
#'
#' @title Unsupervised Maximum Entropy Classifier for Record Linkage
#'
#' @author Adam Struzik
#'
#' @export
mec <- function(A,
                B,
                variables,
                comparators = NULL,
                methods = NULL,
                start_params = NULL,
                set_construction = NULL,
                delta = 0.5,
                eps = 0.05,
                controls_nleqslv = list(),
                controls_kliep = control_kliep(),
                true_matches = NULL) {

  if (!is.null(methods)) {
    stopifnot("`methods` should be a list." =
                is.list(methods))
  }

  if (!is.null(start_params)) {
    stopifnot("`start_params` should be a list." =
                is.list(start_params))
  }

  if (is.null(set_construction)) {
    set_construction <- "size"
  }

  data.table::setDT(A)
  data.table::setDT(B)
  data.table::set(A, j = "a", value = seq_len(nrow(A)))
  data.table::set(B, j = "b", value = seq_len(nrow(B)))
  M <- merge(A, B, by = variables)

  stopifnot("There are no records with perfect agreement on the key variables.
            Please provide relevant datasets." =
              NROW(M) > 0)

  unique_values <- lapply(variables, function(col) {
    length(unique(c(A[[col]], B[[col]])))
  })
  names(unique_values) <- variables

  for (var in variables) {
    if (unique_values[[var]] == 1) {
      data.table::set(A, j = substitute(var), value = NULL)
      data.table::set(B, j = substitute(var), value = NULL)
      variables <- variables[variables != var]
      comparators[[var]] <- NULL
      methods[[var]] <- NULL
      warning(paste("The variable", var, "has only one unique value and has been removed."))
    }
  }

  vectors <- comparison_vectors(A = A,
                                B = B,
                                variables = variables,
                                comparators = comparators)
  Omega <- vectors$Omega
  comparators <- vectors$comparators

  missing_variables <- variables[!(variables %in% names(methods))]
  methods[missing_variables] <- "binary"

  b_vars <- NULL
  cpar_vars <- NULL
  cnonpar_vars <- NULL

  if (any(methods == "binary")) {
    b_vars <- paste0("gamma_", names(which(methods == "binary")))
  }

  if (any(methods == "continuous_parametric")) {
    cpar_vars <- paste0("gamma_", names(which(methods == "continuous_parametric")))
  }

  if (any(methods == "continuous_nonparametric")) {
    cnonpar_vars <- paste0("gamma_", names(which(methods == "continuous_nonparametric")))
  }

  if (is.null(start_params)) {

    start_params <- list()

    if (length(b_vars) > 0) {
      start_params[["binary"]] <- data.table(
        variable = b_vars,
        theta = runif(length(b_vars), min = 0.9)
      )

    }

    if (length(cpar_vars) > 0) {
      start_params[["continuous_parametric"]] <- data.table(
        variable = cpar_vars,
        p_0_M = runif(length(cpar_vars), min = 0.9),
        alpha_M = runif(length(cpar_vars), max = 1),
        beta_M = runif(length(cpar_vars), max = 10)
      )
    }

  }

  ratio_kliep <- NULL

  M <- merge(M, Omega, by = c("a", "b"), all = FALSE)
  M <- M[, colnames(Omega), with = FALSE]
  U <- data.table::fsetdiff(Omega, M)
  n <- NROW(Omega)
  n_M <- NROW(M)
  data.table::set(Omega, j = "ratio", value = 1)

  if (length(b_vars) > 0) {

    b_params <- start_params$binary
    Omega_b <- Omega[, b_vars, with = FALSE]
    M_b <- M[, b_vars, with = FALSE]
    U_b <- U[, b_vars, with = FALSE]

    eta_b <- binary_formula(Omega_b)
    b_params$eta <- eta_b

    b_numerator_list <- lapply(b_vars,
                                    function(col) {
                                      stats::dbinom(x = Omega_b[[col]],
                                                    size = 1,
                                                    prob = as.numeric(b_params[b_params[["variable"]] == col, "theta"]))
                                    })
    b_numerator <- Reduce(`*`, b_numerator_list)
    b_denominator_list <- lapply(b_vars,
                                      function(col) {
                                        stats::dbinom(x = Omega_b[[col]],
                                                      size = 1,
                                                      prob = as.numeric(b_params[b_params[["variable"]] == col, "eta"]))
                                      })
    b_denominator <- Reduce(`*`, b_denominator_list)
    data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * b_numerator / b_denominator)
    data.table::set(Omega, j = "b_denominator", value = b_denominator)

  }

  if (length(cpar_vars) > 0) {

    cpar_params <- start_params$continuous_parametric
    Omega_cpar <- Omega[, cpar_vars, with = FALSE]
    M_cpar <- M[, cpar_vars, with = FALSE]
    U_cpar <- U[, cpar_vars, with = FALSE]

    p_0_U <- p_0_formula(Omega_cpar)
    gamma_plus_U <- gamma_plus_formula(Omega_cpar)
    modified_nleqslv <- purrr::partial(nleqslv::nleqslv, control = controls_nleqslv)
    alpha_U <- alpha_formula(Omega_cpar, modified_nleqslv)
    beta_U <- alpha_U / gamma_plus_U
    cpar_params$p_0_U <- p_0_U
    cpar_params$alpha_U <- alpha_U
    cpar_params$beta_U <- beta_U

    cpar_numerator_list <- lapply(cpar_vars,
                                                   function(col) {
                                                     hurdle_gamma_density(x = Omega_cpar[[col]],
                                                                          p_0 = as.numeric(cpar_params[cpar_params[["variable"]] == col, "p_0_M"]),
                                                                          alpha = as.numeric(cpar_params[cpar_params[["variable"]] == col, "alpha_M"]),
                                                                          beta = as.numeric(cpar_params[cpar_params[["variable"]] == col, "beta_M"]))
                                                   })
    cpar_numerator <- Reduce(`*`, cpar_numerator_list)
    cpar_denominator_list <- lapply(cpar_vars,
                                                     function(col) {
                                                       hurdle_gamma_density(x = Omega_cpar[[col]],
                                                                            p_0 = as.numeric(cpar_params[cpar_params[["variable"]] == col, "p_0_U"]),
                                                                            alpha = as.numeric(cpar_params[cpar_params[["variable"]] == col, "alpha_U"]),
                                                                            beta = as.numeric(cpar_params[cpar_params[["variable"]] == col, "beta_U"]))
                                                     })
    cpar_denominator <- Reduce(`*`, cpar_denominator_list)
    data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * cpar_numerator / cpar_denominator)
    data.table::set(Omega, j = "cpar_denominator", value = cpar_denominator)

  }

  if (length(cnonpar_vars) > 0) {

    Omega_cnonpar <- Omega[, cnonpar_vars, with = FALSE]
    M_cnonpar <- M[, cnonpar_vars, with = FALSE]
    U_cnonpar <- U[, cnonpar_vars, with = FALSE]

    Omega_indexes <- paste0(Omega[["a"]], "_", Omega[["b"]])
    M_indexes <- paste0(M[["a"]], "_", M[["b"]])
    ratio_temp <- as.numeric(Omega$ratio)
    ratio_temp[which(Omega_indexes %in% M_indexes)] <- (Omega$ratio)[which(Omega_indexes %in% M_indexes)] * stats::runif(length(which(Omega_indexes %in% M_indexes)),
                                                                       min = 500, max = 1000)
    ratio_temp[setdiff(1:n, which(Omega_indexes %in% M_indexes))] <- (Omega$ratio)[setdiff(1:n, which(Omega_indexes %in% M_indexes))] * stats::runif(n - length(which(Omega_indexes %in%M_indexes)),
                                                                                                                                                     min = 0.9, max = 1)
    data.table::set(Omega, j = "ratio", value = ratio_temp)

  }

  iter <- 1

  repeat {

    g_est <- pmin(NROW(M) * Omega$ratio / (NROW(M) * (Omega$ratio - 1) + n), 1)
    n_M_old <- n_M
    n_M <- sum(g_est)

    Omega <- Omega[order(-get("ratio")), ]
    M <- data.table("a" = numeric(), "b" = numeric())
    for (var in colnames(Omega)) {
      M[[var]] <- numeric()
    }
    M[["ratio"]] <- numeric()
    used_a <- c()
    used_b <- c()

    for (i in 1:NROW(Omega)) {

      current_a <- Omega$a[i]
      current_b <- Omega$b[i]
      if (!(current_a %in% used_a) && !(current_b %in% used_b)) {
        # M <- rbind(M, Omega[i, c("a", "b", "ratio")])
        M <- rbind(M, Omega[i, ])
        used_a <- c(used_a, current_a)
        used_b <- c(used_b, current_b)
      }

    }

    M <- head(M, round(n_M))
    U <- data.table::fsetdiff(Omega, M)

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

      if (length(cnonpar_vars) == 0) {

        if ((abs(n_M_old - n_M) < delta) || norm(old_params - params, type = "2") < eps) {
          break
        }

      } else {

        if ((abs(n_M_old - n_M) < delta)) {
          break
        }

      }

    }

    data.table::set(Omega, j = "ratio", value = 1)

    if (length(b_vars) > 0) {

      Omega_b <- Omega[, b_vars, with = FALSE]
      M_b <- M[, b_vars, with = FALSE]
      U_b <- U[, b_vars, with = FALSE]

      theta_b_old <- b_params$theta
      theta_b <- binary_formula(M_b)
      b_params$theta <- theta_b

      b_numerator_list <- lapply(b_vars,
                                      function(col) {
                                        stats::dbinom(x = Omega_b[[col]],
                                                      size = 1,
                                                      prob = as.numeric(b_params[b_params[["variable"]] == col, "theta"]))
                                      })
      b_numerator <- Reduce(`*`, b_numerator_list)
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * b_numerator / Omega[["b_denominator"]])
    }

    if (length(cpar_vars) > 0) {

      Omega_cpar <- Omega[, cpar_vars, with = FALSE]
      M_cpar <- M[, cpar_vars, with = FALSE]
      U_cpar <- U[, cpar_vars, with = FALSE]
      p_0_M_old <- cpar_params$p_0_M
      alpha_M_old <- cpar_params$alpha_M
      beta_M_old <- cpar_params$beta_M
      p_0_M <- p_0_formula(M_cpar)
      gamma_plus_M <- gamma_plus_formula(M_cpar)
      alpha_M <- alpha_formula_iterative(M_cpar, modified_nleqslv, beta_M_old)
      beta_M <- alpha_M / gamma_plus_M
      beta_M[is.nan(beta_M)] <- beta_M_old[is.nan(beta_M)]
      cpar_params$p_0_M <- p_0_M
      cpar_params$alpha_M <- alpha_M
      cpar_params$beta_M <- beta_M
      cpar_numerator_list <- lapply(cpar_vars,
                                                     function(col) {
                                                       hurdle_gamma_density(x = Omega_cpar[[col]],
                                                                            p_0 = as.numeric(cpar_params[cpar_params[["variable"]] == col, "p_0_M"]),
                                                                            alpha = as.numeric(cpar_params[cpar_params[["variable"]] == col, "alpha_M"]),
                                                                            beta = as.numeric(cpar_params[cpar_params[["variable"]] == col, "beta_M"]))
                                                     })

      cpar_numerator <- Reduce(`*`, cpar_numerator_list)
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * cpar_numerator / Omega[["cpar_denominator"]])

    }

    if (length(cnonpar_vars) > 0) {

      Omega_cnonpar <- Omega[, cnonpar_vars, with = FALSE]
      M_cnonpar <- M[, cnonpar_vars, with = FALSE]
      U_cnonpar <- U[, cnonpar_vars, with = FALSE]

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

    iter <- iter + 1

  }

  if (!is.null(true_matches)) {

    data.table::setDT(true_matches)

    eval <- evaluation(M, true_matches, n)
    eval_metrics <- unlist(get_metrics(
      TP = eval$TP,
      FP = eval$FP,
      FN = eval$FN,
      TN = eval$TN
    ))
    confusion <- get_confusion(
      TP = eval$TP,
      FP = eval$FP,
      FN = eval$FN,
      TN = eval$TN
    )

  }

  structure(
    list(
      M_est = M[, c("a", "b", "ratio")],
      n_M_est = n_M,
      b_vars = if (length(b_vars) == 0) NULL else b_vars,
      cpar_vars = if (length(cpar_vars) == 0) NULL else cpar_vars,
      cnonpar_vars = if (length(cnonpar_vars) == 0) NULL else cnonpar_vars,
      b_params = if (length(b_vars) == 0) NULL else b_params,
      cpar_params = if (length(cpar_vars) == 0) NULL else cpar_params,
      ratio_kliep = if (is.null(ratio_kliep)) NULL else ratio_kliep,
      variables = variables,
      set_construction = set_construction,
      eval_metrics = if (is.null(true_matches)) NULL else eval_metrics,
      confusion = if (is.null(true_matches)) NULL else confusion
    ),
    class = "mec_rec_lin"
  )

}
