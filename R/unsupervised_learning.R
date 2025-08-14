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
                set_construction = c("size", "flr"),
                delta = 0.5,
                eps = 0.05, # to think
                controls_nleqslv = list(),
                controls_kliep = control_kliep()) {

  if (!is.null(methods)) {
    stopifnot("`methods` should be a list." =
                is.list(methods))
  }

  if (!is.null(start_params)) {
    stopifnot("`start_params` should be a list." =
                is.list(start_params))
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

  binary_variables <- NULL
  continuous_parametric_variables <- NULL
  continuous_nonparametric_variables <- NULL

  if (any(methods == "binary")) {
    binary_variables <- paste0("gamma_", names(which(methods == "binary")))
  }

  if (any(methods == "continuous_parametric")) {
    continuous_parametric_variables <- paste0("gamma_", names(which(methods == "continuous_parametric")))
  }

  if (any(methods == "continuous_nonparametric")) {
    continuous_nonparametric_variables <- paste0("gamma_", names(which(methods == "continuous_nonparametric")))
  }

  if (is.null(start_params)) {

    start_params <- list()

    if (length(binary_variables) > 0) {
      start_params[["binary"]] <- data.table(
        variable = binary_variables,
        theta = runif(length(binary_variables), min = 0.5)#,
        # eta = runif(length(binary_variables))
      )

    }

    if (length(continuous_parametric_variables) > 0) {
      start_params[["continuous_parametric"]] <- data.table(
        variable = continuous_parametric_variables,
        p_0_M = runif(length(continuous_parametric_variables), min = 0.5),
        # p_0_U = runif(length(continuous_parametric_variables)),
        alpha_M = runif(length(continuous_parametric_variables), max = 1),
        # alpha_U = runif(length(continuous_parametric_variables), max = 1000),
        beta_M = runif(length(continuous_parametric_variables), max = 10)#,
        # beta_U = runif(length(continuous_parametric_variables), max = 1000)
      )
    }

  }

  M <- merge(M, Omega, by = c("a", "b"), all = FALSE)
  M <- M[, colnames(Omega), with = FALSE]
  U <- data.table::fsetdiff(Omega, M)
  n <- NROW(Omega)
  n_M <- NROW(M)
  data.table::set(Omega, j = "ratio", value = 1)

  if (length(binary_variables) > 0) {

    binary_params <- start_params$binary
    Omega_binary <- Omega[, binary_variables, with = FALSE]
    M_binary <- M[, binary_variables, with = FALSE]
    U_binary <- U[, binary_variables, with = FALSE]

    eta_binary <- binary_formula(Omega_binary)
    binary_params$eta <- eta_binary

    binary_numerator_list <- lapply(binary_variables,
                                    function(col) {
                                      stats::dbinom(x = Omega_binary[[col]],
                                                    size = 1,
                                                    prob = as.numeric(binary_params[binary_params[["variable"]] == col, "theta"]))
                                    })
    binary_numerator <- Reduce(`*`, binary_numerator_list)
    binary_denominator_list <- lapply(binary_variables,
                                      function(col) {
                                        stats::dbinom(x = Omega_binary[[col]],
                                                      size = 1,
                                                      prob = as.numeric(binary_params[binary_params[["variable"]] == col, "eta"]))
                                      })
    binary_denominator <- Reduce(`*`, binary_denominator_list)
    data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * binary_numerator / binary_denominator)
    data.table::set(Omega, j = "binary_denominator", value = binary_denominator)

  }

  if (length(continuous_parametric_variables) > 0) {

    continuous_parametric_params <- start_params$continuous_parametric
    Omega_continuous_parametric <- Omega[, continuous_parametric_variables, with = FALSE]
    M_continuous_parametric <- M[, continuous_parametric_variables, with = FALSE]
    U_continuous_parametric <- U[, continuous_parametric_variables, with = FALSE]

    p_0_U <- p_0_formula(Omega_continuous_parametric)
    gamma_plus_U <- gamma_plus_formula(Omega_continuous_parametric)
    modified_nleqslv <- purrr::partial(nleqslv::nleqslv, control = controls_nleqslv)
    alpha_U <- alpha_formula(Omega_continuous_parametric, modified_nleqslv)
    beta_U <- alpha_U / gamma_plus_U
    continuous_parametric_params$p_0_U <- p_0_U
    continuous_parametric_params$alpha_U <- alpha_U
    continuous_parametric_params$beta_U <- beta_U

    continuous_parametric_numerator_list <- lapply(continuous_parametric_variables,
                                                   function(col) {
                                                     hurdle_gamma_density(x = Omega_continuous_parametric[[col]],
                                                                          p_0 = as.numeric(continuous_parametric_params[continuous_parametric_params[["variable"]] == col, "p_0_M"]),
                                                                          alpha = as.numeric(continuous_parametric_params[continuous_parametric_params[["variable"]] == col, "alpha_M"]),
                                                                          beta = as.numeric(continuous_parametric_params[continuous_parametric_params[["variable"]] == col, "beta_M"]))
                                                   })
    continuous_parametric_numerator <- Reduce(`*`, continuous_parametric_numerator_list)
    continuous_parametric_denominator_list <- lapply(continuous_parametric_variables,
                                                     function(col) {
                                                       hurdle_gamma_density(x = Omega_continuous_parametric[[col]],
                                                                            p_0 = as.numeric(continuous_parametric_params[continuous_parametric_params[["variable"]] == col, "p_0_U"]),
                                                                            alpha = as.numeric(continuous_parametric_params[continuous_parametric_params[["variable"]] == col, "alpha_U"]),
                                                                            beta = as.numeric(continuous_parametric_params[continuous_parametric_params[["variable"]] == col, "beta_U"]))
                                                     })
    continuous_parametric_denominator <- Reduce(`*`, continuous_parametric_denominator_list)
    data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * continuous_parametric_numerator / continuous_parametric_denominator)
    data.table::set(Omega, j = "continuous_parametric_denominator", value = continuous_parametric_denominator)

  }

  if (length(continuous_nonparametric_variables) > 0) {

    Omega_continuous_nonparametric <- Omega[, continuous_nonparametric_variables, with = FALSE]
    M_continuous_nonparametric <- M[, continuous_nonparametric_variables, with = FALSE]
    U_continuous_nonparametric <- U[, continuous_nonparametric_variables, with = FALSE]

    Omega_indexes <- paste0(Omega[["a"]], "_", Omega[["b"]])
    M_indexes <- paste0(M[["a"]], "_", M[["b"]])
    ratio_temp <- as.numeric(Omega$ratio)
    ratio_temp[which(Omega_indexes %in% M_indexes)] <- (Omega$ratio)[which(Omega_indexes %in% M_indexes)] * stats::runif(length(which(Omega_indexes %in% M_indexes)),
                                                                       min = 1, max = 1.1)
    ratio_temp[setdiff(1:n, which(Omega_indexes %in% M_indexes))] <- (Omega$ratio)[setdiff(1:n, which(Omega_indexes %in% M_indexes))] * stats::runif(n - length(which(Omega_indexes %in%M_indexes)),
                                                                                                                                                     min = 0.9, max = 1)
    data.table::set(Omega, j = "ratio", value = ratio_temp)

    # ratio_kliep <- do.call(
    #   densityratio::kliep,
    #   c(list(
    #     df_numerator = M_continuous_nonparametric,
    #     df_denominator = U_continuous_nonparametric
    #   ),
    #   controls_kliep)
    # )
    #
    # data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * stats::predict(ratio_kliep, Omega_continuous_nonparametric))

  }

  repeat {

    g_est <- pmin(n_M * Omega$ratio / (n_M * (Omega$ratio - 1) + n), 1)
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
    data.table::set(Omega, j = "ratio", value = 1)

    if (length(binary_variables) > 0) {

      M_binary <- M[, binary_variables, with = FALSE]
      U_binary <- U[, binary_variables, with = FALSE]

      theta_binary_old <- binary_params$theta
      theta_binary <- binary_formula(M_binary)
      binary_params$theta <- theta_binary

      binary_numerator_list <- lapply(binary_variables,
                                      function(col) {
                                        stats::dbinom(x = Omega_binary[[col]],
                                                      size = 1,
                                                      prob = as.numeric(binary_params[binary_params[["variable"]] == col, "theta"]))
                                      })
      binary_numerator <- Reduce(`*`, binary_numerator_list)
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * binary_numerator / Omega[["binary_denominator"]])
    }

    if (length(continuous_parametric_variables) > 0) {

      M_continuous_parametric <- M[, continuous_parametric_variables, with = FALSE]
      U_continuous_parametric <- U[, continuous_parametric_variables, with = FALSE]
      p_0_M_old <- continuous_parametric_params$p_0_M
      alpha_M_old <- continuous_parametric_params$alpha_M
      beta_M_old <- continuous_parametric_params$beta_M
      p_0_M <- p_0_formula(M_continuous_parametric)
      gamma_plus_M <- gamma_plus_formula(M_continuous_parametric)
      alpha_M <- alpha_formula_iterative(M_continuous_parametric, modified_nleqslv, beta_M_old)
      beta_M <- alpha_M / gamma_plus_M
      beta_M[is.nan(beta_M)] <- beta_M_old[is.nan(beta_M)]
      continuous_parametric_params$p_0_M <- p_0_M
      continuous_parametric_params$alpha_M <- alpha_M
      continuous_parametric_params$beta_M <- beta_M
      continuous_parametric_numerator_list <- lapply(continuous_parametric_variables,
                                                     function(col) {
                                                       hurdle_gamma_density(x = Omega_continuous_parametric[[col]],
                                                                            p_0 = as.numeric(continuous_parametric_params[continuous_parametric_params[["variable"]] == col, "p_0_M"]),
                                                                            alpha = as.numeric(continuous_parametric_params[continuous_parametric_params[["variable"]] == col, "alpha_M"]),
                                                                            beta = as.numeric(continuous_parametric_params[continuous_parametric_params[["variable"]] == col, "beta_M"]))
                                                     })
      continuous_parametric_numerator <- Reduce(`*`, continuous_parametric_numerator_list)
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * continuous_parametric_numerator / Omega[["continuous_parametric_denominator"]])

    }

    if (length(continuous_nonparametric_variables) > 0) {

      M_continuous_nonparametric <- M[, continuous_nonparametric_variables, with = FALSE]
      U_continuous_nonparametric <- U[, continuous_nonparametric_variables, with = FALSE]

      ratio_kliep <- do.call(
        densityratio::kliep,
        c(list(
          df_numerator = M_continuous_nonparametric,
          df_denominator = U_continuous_nonparametric
        ),
        controls_kliep)
      )

      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * stats::predict(ratio_kliep, Omega_continuous_nonparametric))

    }

  }

}
