#' @import data.table
#' @importFrom nleqslv nleqslv
#' @importFrom purrr partial
#' @importFrom densityratio kliep
#'
#' @title Train a Record Linkage Model
#'
#' @author Adam Struzik
#'
#' @description Trains a supervised record linkage model using probability or density ratio estimation.
#'
#' @param A A duplicate-free `data.frame` or `data.table`.
#' @param B A duplicate-free `data.frame` or `data.table`.
#' @param matches A `data.frame` or `data.table` indicating known matches.
#' @param variables A character vector of key variables used to create comparison vectors.
#' @param comparators A named list of functions for comparing pairs of records.
#' @param methods A named list of methods used for estimation (`"binary"`, `"continuous_parametric"` or `"continuous_nonparametric"`).
#' @param controls_nleqslv Controls passed to the \link[nleqslv]{nleqslv} function (only if the `"continuous_parametric"` method has been chosen for at least one variable).
#' @param controls_kliep Controls passed to the \link[densityratio]{kliep} function (only if the `"continuous_nonparametric"` method has been chosen for at least one variable).
#'
#' @return
#' Returns a list containing:\cr
#' \itemize{
#' \item{`binary_variables` -- a character vector of variables chosen for the `"binary"` method (with the prefix `"gamma_"`),}
#' \item{`continuous_parametric_variables` -- a character vector of variables chosen for the `"continuous_parametric"` method (with the prefix `"gamma_"`),}
#' \item{`continuous_nonparametric_variables` -- a character vector of variables chosen for the `"continuous_nonparametric"` method (with the prefix `"gamma_"`),}
#' \item{`binary_params` -- parameters estimated using the `"binary"` method,}
#' \item{`continuous_parametric_params` -- parameters estimated using the `"continuous_parametric"` method,}
#' \item{`ratio_kliep` -- a result of the \link[densityratio]{kliep} function,}
#' \item{`ml_model` -- here `NULL`,}
#' \item{`pi_est` -- a prior probability of matching,}
#' \item{`variables` -- a character vector of key variables used for comparison,}
#' \item{`comparators` -- a list of functions used to compare pairs of records,}
#' \item{`methods` -- a list of methods used for estimation.}
#' }
#'
#' @examples
#' df_1 <- data.frame(
#'   "name" = c("John", "Emily", "Mark", "Anna", "David"),
#'   "surname" = c("Smith", "Johnson", "Taylor", "Williams", "Brown")
#' )
#' df_2 <- data.frame(
#'   "name" = c("Jon", "Emely", "Marc", "Michael"),
#'   "surname" = c("Smitth", "Jonson", "Tailor", "Henderson")
#' )
#' comparators <- list("name" = reclin2::cmp_jarowinkler(),
#'                     "surname" = reclin2::cmp_jarowinkler())
#' matches <- data.frame("a" = 1:3, "b" = 1:3)
#' methods <- list("name" = "continuous_nonparametric",
#'                 "surname" = "continuous_nonparametric")
#' model <- train_rec_lin(A = df_1, B = df_2, matches = matches,
#'                        variables = c("name", "surname"),
#'                        comparators = comparators,
#'                        methods = methods,
#'                        controls_kliep = control_kliep(nfold = 3))
#' model
#'
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

#' @title Create a Custom Record Linkage Model
#'
#' @author Adam Struzik
#'
#' @description
#' Creates a supervised record linkage model using a custom machine learning (ML) classifier.
#'
#' @param ml_model A trained ML model that predicts the probability of a match based on comparison vectors.
#' @param vectors An object of class `comparison_vectors` (a result of the `comparison_vectors` function), used for training the `ml_model`.
#'
#' @return
#' Returns a list containing:\cr
#' \itemize{
#' \item{`binary_variables` -- here `NULL`,}
#' \item{`continuous_parametric_variables` -- here `NULL`,}
#' \item{`continuous_nonparametric_variables` -- here `NULL`,}
#' \item{`binary_params` -- here `NULL`,}
#' \item{`continuous_parametric_params` -- here `NULL`,}
#' \item{`ratio_kliep` -- here `NULL`,}
#' \item{`ml_model` -- ML model used for creating the record linkage model,}
#' \item{`pi_est` -- a prior probability of matching,}
#' \item{`variables` -- a character vector of key variables used for comparison,}
#' \item{`comparators` -- a list of functions used to compare pairs of records,}
#' \item{`methods` -- here `NULL`.}
#' }
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
      methods = NULL
    ),
    class = "rec_lin_model"
  )

}
