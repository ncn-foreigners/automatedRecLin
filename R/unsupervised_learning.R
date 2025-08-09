#' @import data.table
#' @importFrom stats runif
#' @importFrom purrr partial
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
        theta = runif(length(binary_variables))#,
        # eta = runif(length(binary_variables))
      )

    }

    if (length(continuous_parametric_variables) > 0) {
      start_params[["continuous_parametric"]] <- data.table(
        variable = continuous_parametric_variables,
        p_0_M = runif(length(continuous_parametric_variables)),
        # p_0_U = runif(length(continuous_parametric_variables)),
        alpha_M = runif(length(continuous_parametric_variables), max = 100),
        # alpha_U = runif(length(continuous_parametric_variables), max = 1000),
        beta_M = runif(length(continuous_parametric_variables), max = 100)#,
        # beta_U = runif(length(continuous_parametric_variables), max = 1000)
      )
    }

  }

  M <- merge(M, Omega, by = c("a", "b"), all = FALSE)
  M <- M[, colnames(Omega), with = FALSE]
  U <- data.table::fsetdiff(Omega, M)
  n <- NROW(Omega)
  n_M <- NROW(M)

  if (length(binary_variables) > 0) {

    binary_params <- start_params$binary
    Omega_binary <- Omega[, binary_variables, with = FALSE]
    M_binary <- M[, binary_variables, with = FALSE]
    U_binary <- U[, binary_variables, with = FALSE]

    eta_binary <- binary_formula(Omega_binary)
    binary_params$eta <- eta_binary

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

  }

  if (length(continuous_nonparametric_variables) > 0) {

    # Omega_continuous_nonparametric <- Omega[, continuous_nonparametric_variables, with = FALSE]
    M_continuous_nonparametric <- M[, continuous_nonparametric_variables, with = FALSE]
    U_continuous_nonparametric <- U[, continuous_nonparametric_variables, with = FALSE]

    # ratio_kliep <- do.call(
    #   densityratio::kliep,
    #   c(list(
    #     df_numerator = M_continuous_nonparametric,
    #     df_denominator = U_continuous_nonparametric
    #   ),
    #   controls_kliep)
    # )

  }

  repeat {

    if (length(binary_variables) > 0) {

      theta_binary_old <- binary_params$theta
      theta_binary <- binary_formula(M_binary)
      binary_params$theta <- theta_binary

    }

    if (length(continuous_parametric_variables) > 0) {

      p_0_M_old <- continuous_parametric_params$p_0_M
      alpha_M_old <- continuous_parametric_params$alpha_M
      beta_M_old <- continuous_parametric_params$beta_M
      p_0_M <- p_0_formula(M_continuous_parametric)
      gamma_plus_M <- gamma_plus_formula(M_continuous_parametric)
      alpha_M <- alpha_formula_iterative(M_continuous_parametric, modified_nleqslv, beta_M_old)

    }

    if (length(continuous_nonparametric_variables) > 0) {

      ratio_kliep <- do.call(
        densityratio::kliep,
        c(list(
          df_numerator = M_continuous_nonparametric,
          df_denominator = U_continuous_nonparametric
        ),
        controls_kliep)
      )

    }

  }

}
