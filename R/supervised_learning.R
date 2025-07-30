#' @import data.table
#' @importFrom nleqslv nleqslv
#' @importFrom purrr partial
#' @importFrom densityratio kliep
#'
#' @title Train a record linkage model
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

  stopifnot("`methods` should be a list." =
            is.list(methods))

  missing_variables <- variables[!(variables %in% names(methods))]
  methods[missing_variables] <- "binary"

  M <- Omega[match == 1, ]
  U <- Omega[match == 0, ]

  n <- NROW(Omega)
  n_M <- NROW(M)
  pi_est <- n_M / n

  binary_variables <- NULL
  continuous_parametric_variables <- NULL
  continuous_nonparametric_variables <- NULL

  binary_params <- NULL
  continuous_parametric_params <- NULL
  ratio_kliep <- NULL

  if (any(methods == "binary")) {

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

  if (any(methods == "continuous_parametric")) {

    continuous_parametric_variables <- paste0("gamma_", names(which(methods == "continuous_parametric")))
    M_continuous_parametric <- M[, ..continuous_parametric_variables]
    U_continuous_parametric <- U[, ..continuous_parametric_variables]
    p_0_M <- p_0_formula(M_continuous_parametric)
    p_0_U <- p_0_formula(U_continuous_parametric)
    gamma_plus_M <- gamma_plus_formula(M_continuous_parametric)
    gamma_plus_U <- gamma_plus_formula(U_continuous_parametric)
    modified_nleqslv <- purrr::partial(nleqslv::nleqslv, control = controls_nleqslv)
    alpha_M <- alpha_formula(M_continuous_parametric, modified_nleqslv)
    alpha_U <- alpha_formula(U_continuous_parametric, modified_nleqslv)
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

  if (any(methods == "continuous_nonparametric")) {

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

  structure(
    list(
      binary_variables = if (is.null(binary_variables)) NULL else binary_variables,
      continuous_parametric_variables = if (is.null(continuous_parametric_variables)) NULL else continuous_parametric_variables,
      continuous_nonparametric_variables = if (is.null(continuous_nonparametric_variables)) NULL else continuous_nonparametric_variables,
      binary_params = if (is.null(binary_params)) NULL else binary_params,
      continuous_parametric_params = if (is.null(continuous_parametric_params)) NULL else continuous_parametric_params,
      ratio_kliep = if (is.null(ratio_kliep)) NULL else ratio_kliep,
      ml_model = NULL,
      pi_est = pi_est,
      variables = variables,
      comparators = comparators,
      methods = methods
    ),
    class = "rec_lin_model"
  )

}

#' @title Create a custom record linkage model
#'
#' @export
custom_rec_lin_model <- function(ml_model, vectors) {

  stopifnot("`vectors` must be of class `comparison_vectors`" =
            is(vectors, "comparison_vectors"))

  Omega <- vectors$Omega

  M <- Omega[match == 1, ]
  U <- Omega[match == 0, ]

  n <- NROW(Omega)
  n_M <- NROW(M)
  pi_est <- n_M / n

  variables <- vectors$variables
  K <- length(variables)

  methods <- replicate(K, "custom", simplify = FALSE)
  names(methods) <- variables

  structure(
    list(
      binary_variables = NULL,
      continuous_parametric_variables = NULL,
      continuous_nonparametric_variables = NULL,
      binary_params = NULL,
      continuous_parametric_params = NULL,
      ratio_kliep = NULL,
      ml_model = ml_model,
      pi_est = pi_est,
      variables = variables,
      comparators = vectors$comparators,
      methods = methods
    ),
    class = "rec_lin_model"
  )

}
