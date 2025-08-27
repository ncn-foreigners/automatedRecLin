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
#' \item{`b_vars` -- a character vector of variables used for the `"binary"` method (with the prefix `"gamma_"`),}
#' \item{`cpar_vars` -- a character vector of variables used for the `"continuous_parametric"` method (with the prefix `"gamma_"`),}
#' \item{`cnonpar_vars` -- a character vector of variables used for the `"continuous_nonparametric"` method (with the prefix `"gamma_"`),}
#' \item{`b_params` -- parameters estimated using the `"binary"` method,}
#' \item{`cpar_params` -- parameters estimated using the `"continuous_parametric"` method,}
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
#'   "name" = c("John", "Emely", "Marc", "Michael"),
#'   "surname" = c("Smith", "Jonson", "Tailor", "Henderson")
#' )
#' comparators <- list("name" = jarowinkler_complement(),
#'                     "surname" = jarowinkler_complement())
#' matches <- data.frame("a" = 1:3, "b" = 1:3)
#' methods <- list("name" = "continuous_nonparametric",
#'                 "surname" = "continuous_nonparametric")
#' model <- train_rec_lin(A = df_1, B = df_2, matches = matches,
#'                        variables = c("name", "surname"),
#'                        comparators = comparators,
#'                        methods = methods)
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

  data.table::setDT(A)
  data.table::setDT(B)

  stopifnot("There are no records with perfect agreement on the key variables. Please provide relevant datasets." =
              NROW(merge(A, B, by = variables)) > 0)

  if (!is.null(methods)) {
    stopifnot("`methods` should be a list." =
                is.list(methods))
  } else {
    methods <- list()
  }

  missing_variables <- variables[!(variables %in% names(methods))]
  methods[missing_variables] <- "binary"
  methods <- methods[variables]

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

  if (is.null(prob_ratio)) {
    prob_ratio <- "2"
  }

  if (any(methods %in% c("binary", "continuous_parametric"))) {

    data.table::setDT(matches)

    if (prob_ratio == "2") {
      A_temp <- data.table::copy(A)
      B_temp <- data.table::copy(B)
      data.table::set(A_temp, j = "a", value = seq_len(nrow(A)))
      data.table::set(B_temp, j = "b", value = seq_len(nrow(B)))
      indexes <- data.table::CJ(a = A_temp[["a"]], b = B_temp[["b"]])
      data.table::setkey(indexes, NULL)
      indexes_U <- data.table::fsetdiff(indexes, matches)
      for (var in variables) {
        common_values <- sum((A_temp[[var]])[indexes_U$a] == (B_temp[[var]])[indexes_U$b])
        if (common_values == 0) {
          prob_ratio <- "1"
          warning("Some variables lack common values between the unmatches. Switching the probability ratio to \"1\".")
          break
        }
      }
    }

    # if (prob_ratio == "1") {
    #   common_values <- sapply(variables, function(col) {
    #     length(intersect(A[[col]], B[[col]]))
    #   })
    #
    #   stopifnot("Some variables lack common values in both datasets.
    #             Please provide relevant datasets." =
    #               all(common_values > 0))
    # }

  }

  vectors <- comparison_vectors(A = A,
                                B = B,
                                variables = variables,
                                comparators = comparators,
                                matches = matches)
  Omega <- vectors$Omega
  comparators <- vectors$comparators

  M <- Omega[match == 1, ]
  U <- Omega[match == 0, ]

  n <- NROW(Omega)
  n_M <- NROW(M)
  pi_est <- n_M / n

  b_vars <- NULL
  cpar_vars <- NULL
  cnonpar_vars <- NULL

  b_params <- NULL
  cpar_params <- NULL
  ratio_kliep <- NULL

  if (any(methods == "binary")) {

    b_vars <- paste0("gamma_", names(which(methods == "binary")))
    M_b <- M[, b_vars, with = FALSE]
    theta_b <- binary_formula(M_b)
    if (prob_ratio == "1") {
      Omega_b <- Omega[, b_vars, with = FALSE]
      eta_b <- binary_formula(Omega_b)
    } else if (prob_ratio == "2") {
      U_b <- U[, b_vars, with = FALSE]
      eta_b <- binary_formula(U_b)
    }

    b_params <- data.table(
      variable = b_vars,
      theta = theta_b,
      eta = eta_b
    )

  }

  if (any(methods == "continuous_parametric")) {

    cpar_vars <- paste0("gamma_", names(which(methods == "continuous_parametric")))
    modified_nleqslv <- purrr::partial(nleqslv::nleqslv, control = controls_nleqslv)
    M_cpar <- M[, cpar_vars, with = FALSE]
    p_0_M <- p_0_formula(M_cpar)
    gamma_plus_M <- gamma_plus_formula(M_cpar)
    alpha_M <- alpha_formula(M_cpar, modified_nleqslv)
    beta_M <- alpha_M / gamma_plus_M
    if (prob_ratio == "1") {
      Omega_cpar <- Omega[, cpar_vars, with = FALSE]
      p_0_Omega <- p_0_formula(Omega_cpar)
      gamma_plus_Omega <- gamma_plus_formula(Omega_cpar)
      alpha_Omega <- alpha_formula(Omega_cpar, modified_nleqslv)
      beta_Omega <- alpha_Omega / gamma_plus_Omega

      cpar_params <- data.table(
        variable = cpar_vars,
        p_0_M = p_0_M,
        p_0_Omega = p_0_Omega,
        alpha_M = alpha_M,
        alpha_Omega = alpha_Omega,
        beta_M = beta_M,
        beta_Omega = beta_Omega
      )
    } else if (prob_ratio == "2") {
      U_cpar <- U[, cpar_vars, with = FALSE]
      p_0_U <- p_0_formula(U_cpar)
      gamma_plus_U <- gamma_plus_formula(U_cpar)
      alpha_U <- alpha_formula(U_cpar, modified_nleqslv)
      beta_U <- alpha_U / gamma_plus_U

      cpar_params <- data.table(
        variable = cpar_vars,
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

    cnonpar_vars <- paste0("gamma_", names(which(methods == "continuous_nonparametric")))
    M_cnonpar <- M[, cnonpar_vars, with = FALSE]
    if (prob_ratio == "1") {
      Omega_cnonpar <- Omega[, cnonpar_vars, with = FALSE]

      ratio_kliep <- do.call(
        densityratio::kliep,
        c(list(
          df_numerator = M_cnonpar,
          df_denominator = Omega_cnonpar
        ),
        controls_kliep)
      )
    } else if (prob_ratio == "2") {
      U_cnonpar <- U[, cnonpar_vars, with = FALSE]

      ratio_kliep <- do.call(
        densityratio::kliep,
        c(list(
          df_numerator = M_cnonpar,
          df_denominator = U_cnonpar
        ),
        controls_kliep)
      )
    }

  }

  structure(
    list(
      b_vars = if (is.null(b_vars)) NULL else b_vars,
      cpar_vars = if (is.null(cpar_vars)) NULL else cpar_vars,
      cnonpar_vars = if (is.null(cnonpar_vars)) NULL else cnonpar_vars,
      b_params = if (is.null(b_params)) NULL else b_params,
      cpar_params = if (is.null(cpar_params)) NULL else cpar_params,
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
#' \item{`b_vars` -- here `NULL`,}
#' \item{`cpar_vars` -- here `NULL`,}
#' \item{`cnonpar_vars` -- here `NULL`,}
#' \item{`b_params` -- here `NULL`,}
#' \item{`cpar_params` -- here `NULL`,}
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
#'   comparators <- list("name" = jarowinkler_complement(),
#'                       "surname" = jarowinkler_complement())
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
      b_vars = NULL,
      cpar_vars = NULL,
      cnonpar_vars = NULL,
      b_params = NULL,
      cpar_params = NULL,
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

