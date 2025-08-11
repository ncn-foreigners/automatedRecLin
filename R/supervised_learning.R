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
#' @param prob_ratio Probability ratio type (`"1"` or `"2"`).
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
#' @references
#' Lee, D., Zhang, L.-C. and Kim, J. K. (2022). Maximum entropy classification for record linkage.
#' Survey Methodology, Statistics Canada, Catalogue No. 12-001-X, Vol. 48, No. 1.
#'
#' Vo, T. H., Chauvet, G., Happe, A., Oger, E., Paquelet, S., and Garès, V. (2023).
#' Extending the Fellegi-Sunter record linkage model for mixed-type data with application to the French national health data system.
#' Computational Statistics & Data Analysis, 179, 107656.
#'
#' Sugiyama, M., Suzuki, T., Nakajima, S. et al. Direct importance estimation for covariate shift adaptation.
#' Ann Inst Stat Math 60, 699–746 (2008). \doi{10.1007/s10463-008-0197-x}
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
#'                        prob_ratio = "2",
#'                        controls_kliep = control_kliep(nfold = 3))
#' model
#' @export
train_rec_lin <- function(
    A,
    B,
    matches,
    variables,
    comparators = NULL,
    methods = NULL,
    prob_ratio = NULL,
    controls_nleqslv = list(),
    controls_kliep = control_kliep()) {

  if (!is.null(methods)) {
    stopifnot("`methods` should be a list." =
                is.list(methods))
  }

  if (is.null(prob_ratio)) {
    prob_ratio <- "2"
  }

  vectors <- comparison_vectors(A = A,
                                B = B,
                                variables = variables,
                                comparators = comparators,
                                matches = matches)
  Omega <- vectors$Omega
  comparators <- vectors$comparators

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
    # M_binary <- M[, ..binary_variables]
    # U_binary <- U[, ..binary_variables]
    M_binary <- M[, binary_variables, with = FALSE]
    theta_binary <- binary_formula(M_binary)
    if (prob_ratio == "1") {
      Omega_binary <- Omega[, binary_variables, with = FALSE]
      eta_binary <- binary_formula(Omega_binary)
    } else if (prob_ratio == "2") {
      U_binary <- U[, binary_variables, with = FALSE]
      eta_binary <- binary_formula(U_binary)
    }

    binary_params <- data.table(
      variable = binary_variables,
      theta = theta_binary,
      eta = eta_binary
    )

  }

  if (any(methods == "continuous_parametric")) {

    continuous_parametric_variables <- paste0("gamma_", names(which(methods == "continuous_parametric")))
    # M_continuous_parametric <- M[, ..continuous_parametric_variables]
    # U_continuous_parametric <- U[, ..continuous_parametric_variables]
    modified_nleqslv <- purrr::partial(nleqslv::nleqslv, control = controls_nleqslv)
    M_continuous_parametric <- M[, continuous_parametric_variables, with = FALSE]
    p_0_M <- p_0_formula(M_continuous_parametric)
    gamma_plus_M <- gamma_plus_formula(M_continuous_parametric)
    alpha_M <- alpha_formula(M_continuous_parametric, modified_nleqslv)
    beta_M <- alpha_M / gamma_plus_M
    if (prob_ratio == "1") {
      Omega_continuous_parametric <- Omega[, continuous_parametric_variables, with = FALSE]
      p_0_Omega <- p_0_formula(Omega_continuous_parametric)
      gamma_plus_Omega <- gamma_plus_formula(Omega_continuous_parametric)
      alpha_Omega <- alpha_formula(Omega_continuous_parametric, modified_nleqslv)
      beta_Omega <- alpha_Omega / gamma_plus_Omega

      continuous_parametric_params <- data.table(
        variable = continuous_parametric_variables,
        p_0_M = p_0_M,
        p_0_Omega = p_0_Omega,
        alpha_M = alpha_M,
        alpha_Omega = alpha_Omega,
        beta_M = beta_M,
        beta_Omega = beta_Omega
      )
    } else if (prob_ratio == "2") {
      U_continuous_parametric <- U[, continuous_parametric_variables, with = FALSE]
      p_0_U <- p_0_formula(U_continuous_parametric)
      gamma_plus_U <- gamma_plus_formula(U_continuous_parametric)
      alpha_U <- alpha_formula(U_continuous_parametric, modified_nleqslv)
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

  }

  if (any(methods == "continuous_nonparametric")) {

    continuous_nonparametric_variables <- paste0("gamma_", names(which(methods == "continuous_nonparametric")))
    # M_continuous_nonparametric <- M[, ..continuous_nonparametric_variables]
    # U_continuous_nonparametric <- U[, ..continuous_nonparametric_variables]
    M_continuous_nonparametric <- M[, continuous_nonparametric_variables, with = FALSE]
    if (prob_ratio == "1") {
      Omega_continuous_nonparametric <- Omega[, continuous_nonparametric_variables, with = FALSE]

      ratio_kliep <- do.call(
        densityratio::kliep,
        c(list(
          df_numerator = M_continuous_nonparametric,
          df_denominator = Omega_continuous_nonparametric
        ),
        controls_kliep)
      )
    } else if (prob_ratio == "2") {
      U_continuous_nonparametric <- U[, continuous_nonparametric_variables, with = FALSE]

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

#' @importFrom methods is
#'
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
#' @examples
#' if (requireNamespace("xgboost", quietly = TRUE)) {
#'   df_1 <- data.frame(
#'     "name" = c("John", "Emily", "Mark", "Anna", "David"),
#'     "surname" = c("Smith", "Johnson", "Taylor", "Williams", "Brown")
#'   )
#'   df_2 <- data.frame(
#'     "name" = c("Jon", "Emely", "Marc", "Michael"),
#'     "surname" = c("Smitth", "Jonson", "Tailor", "Henderson")
#'   )
#'   comparators <- list("name" = reclin2::cmp_jarowinkler(),
#'                       "surname" = reclin2::cmp_jarowinkler())
#'   matches <- data.frame("a" = 1:3, "b" = 1:3)
#'   vectors <- comparison_vectors(A = df_1, B = df_2, variables = c("name", "surname"),
#'                                comparators = comparators, matches = matches)
#'   train_data <- xgboost::xgb.DMatrix(
#'     data = as.matrix(vectors$Omega[, c("gamma_name", "gamma_surname")]),
#'     label = vectors$Omega$match
#'   )
#'   params <- list(objective = "binary:logistic",
#'                  eval_metric = "logloss")
#'   model_xgb <- xgboost::xgboost(data = train_data, params = params,
#'                                 nrounds = 50, verbose = 0)
#'   custom_xgb_model <- custom_rec_lin_model(model_xgb, vectors)
#'   custom_xgb_model
#' }
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

