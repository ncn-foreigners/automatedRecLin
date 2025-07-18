#' @import data.table
#' @importFrom nleqslv nleqslv
#' @importFrom purrr partial
#' @importFrom densityratio kliep
#'
#' @title Train a record linkage classifier
#'
#' @export
train_rec_lin <- function(
    A,
    B,
    matches,
    variables,
    comparators = NULL,
    methods = NULL,
    controls_nleqslv = list(),
    controls_kliep = control_kliep()) {

  vectors <- comparison_vectors(A = A,
                                B = B,
                                variables = variables,
                                comparators = comparators,
                                matches = matches)
  Omega <- vectors$Omega
  comparators <- vectors$comparators

  if(!is.null(methods)) {
    stopifnot("`methods` should be a list." =
              is.list(methods))
  }

  missing_variables <- variables[!(variables %in% names(methods))]
  methods[missing_variables] <- "binary"

  M <- Omega[match == 1, ]
  U <- Omega[match ==0, ]

  n <- NROW(Omega)
  n_M <- NROW(M)
  pi_est <- n_M / n

  if (sum(methods == "binary") > 0) {
    binary_variables <- paste0("gamma_", names(which(methods == "binary")))
    M_binary <- M[, ..binary_variables]
    U_binary <- U[, ..binary_variables]
    theta_binary <- binary_formula(M_binary)
    eta_binary <- binary_formula(U_binary)
    binary_params <- data.table(
      variable = binary_variables,
      theta = theta_binary,
      eta = eta_binary
    )
  }

  if (sum(methods == "continuous_parametric") > 0) {
    continuous_parametric_variables <- paste0("gamma_", names(which(methods == "continuous_parametric")))
    M_continuous_parametric <- M[, ..continuous_parametric_variables]
    U_continuous_parametric <- U[, ..continuous_parametric_variables]
    p_0_M <- p_0_formula(M_continuous_parametric)
    p_0_U <- p_0_formula(U_continuous_parametric)
    gamma_plus_M <- gamma_plus_formula(M_continuous_parametric)
    gamma_plus_U <- gamma_plus_formula(U_continuous_parametric)
    temp_function <- purrr::partial(nleqslv::nleqslv, control = controls_nleqslv)
    alpha_M <- alpha_formula(M_continuous_parametric, temp_function)
    alpha_U <- alpha_formula(U_continuous_parametric, temp_function)
    beta_M <- alpha_M / gamma_plus_M
    beta_U <- alpha_U / gamma_plus_U
    continuous_parametric_params <- data.table(
      variable = continuous_parametric_variables,
      p_0_M = p_0_M,
      p_0_U = p_0_U,
      alpha_M = alpha_M,
      alpha_U = alpha_U,
      beta_M = beta_M,
      beta_U = beta_U
    )
  }

  if (sum(methods == "continuous_nonparametric") > 0) {
    continuous_nonparametric_variables <- paste0("gamma_", names(which(methods == "continuous_nonparametric")))
    M_continuous_nonparametric <- M[, ..continuous_nonparametric_variables]
    U_continuous_nonparametric <- U[, ..continuous_nonparametric_variables]
    ratio_kliep <- do.call(
      densityratio::kliep,
      c(list(
        df_numerator = M_continuous_nonparametric,
        df_denominator = U_continuous_nonparametric
      ),
      controls_kliep)
    )
  }

  return(ratio_kliep)

}
