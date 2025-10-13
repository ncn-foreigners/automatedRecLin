#' @import data.table
#' @importFrom nleqslv nleqslv
#' @importFrom purrr partial
#' @importFrom densityratio kliep
#'
#' @title Train a Record Linkage Model
#'
#' @author Adam Struzik
#'
#' @description Trains a supervised record linkage model using probability or density ratio estimation,
#' based on [Lee et al. (2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm),
#' with several extensions.
#'
#' @param A A duplicate-free `data.frame` or `data.table`.
#' @param B A duplicate-free `data.frame` or `data.table`.
#' @param matches A `data.frame` or `data.table` indicating known matches.
#' @param variables A character vector of key variables used to create comparison vectors.
#' @param comparators A named list of functions for comparing pairs of records.
#' @param methods A named list of methods used for estimation (`"binary"`, `"continuous_parametric"` or `"continuous_nonparametric"`).
#' @param prob_ratio Probability ratio type (`"1"` or `"2"`).
#' @param nonpar_hurdle Logical indicating whether to use a hurdle model or not
#' (used only if the `"continuous_nonparametric"` method has been chosen for at least one variable).
#' @param controls_nleqslv Controls passed to the \link[nleqslv]{nleqslv} function (only if the `"continuous_parametric"` method has been chosen for at least one variable).
#' @param controls_kliep Controls passed to the \link[densityratio]{kliep} function (only if the `"continuous_nonparametric"` method has been chosen for at least one variable).
#'
#' @details
#' Consider two datasets: \eqn{A} and \eqn{B}.
#' Let the bipartite comparison space \eqn{\Omega = A \times B} consist of
#' matches \eqn{M} and non-matches \eqn{U} between the records in files
#' \eqn{A} and \eqn{B}. For any pair of records \eqn{(a,b) \in \Omega},
#' let \eqn{\pmb{\gamma}_{ab} = (\gamma_{ab}^1,\gamma_{ab}^2,
#' \ldots,\gamma_{ab}^K)'} be the comparison vector between
#' a set of key variables. The original MEC algorithm uses the binary
#' comparison function to evaluate record pairs across two datasets.
#' However, this approach may be insufficient when handling datasets
#' with frequent errors across multiple variables.
#'
#' We propose the use of continuous comparison functions to address
#' the limitations of binary comparison methods. We consider every
#' semi-metric, i.e., a function \eqn{d: A \times B \to \mathbb{R}},
#' satisfying the following conditions:\cr
#' \enumerate{
#' \item{\eqn{d(x,y) \geq 0},}
#' \item{\eqn{d(x,y) = 0} if and only if \eqn{x = y},}
#' \item{\eqn{d(x,y) = d(y,x)}.}
#' }
#' For example, we can use \eqn{1 - \text{Jaro-Winkler distance}} for character variables
#' (which is implemented in the `automatedRecLin` package as the `jarowinkler_complement` function)
#' or the Euclidean distance for numerical variables. The `automatedRecLin` package allows the use of
#' a different comparison function for each key variable (which should be specified
#' as a list in the `comparators` argument). The default function
#' for each key variable is \link[reclin2]{cmp_identical}
#' (the binary comparison function).
#'
#' The `train_rec_lin` function is used to train a record linkage model,
#' when \eqn{M} and \eqn{U} are known (which might later serve as a classifier
#' for pairs outside \eqn{\Omega}). It offers different approaches to estimate the
#' probability/density ratio between matches and non-matches, which should be
#' specified as a list in the methods argument. The method suitable for the binary
#' comparison function is `"binary"`, which is also the default method for each
#' variable.
#'
#' For the continuous semi-metrics we suggest the usage
#' of `"continuous_parametric"` or `"continuous_nonparametric"`
#' method. The `"continuous_parametric"` method assumes that
#' \eqn{\gamma_{ab}^k|M} and \eqn{\gamma_{ab}^k|U} follow
#' hurdle Gamma distributions. The density function of a hurdle
#' Gamma distribution is characterized by three parameters
#' \eqn{p_0 \in (0,1)} and \eqn{\alpha, \beta > 0} as follows:
#' \deqn{
#' f(x;p_0,\alpha,\beta) = p_0^{\mathbb{I}(x = 0)}[(1 - p_0) v(x;\alpha,\beta)]^{\mathbb{I}(x > 0)},
#' }
#' where
#' \deqn{
#' v(x;\alpha,\beta) = \frac{\beta^{\alpha} x^{\alpha - 1} \exp(-\beta x)}
#' {\Gamma(\alpha)}
#' }
#' is the density function of a Gamma distribution
#' (for details see [Vo et al. (2023)](https://ideas.repec.org/a/eee/csdana/v179y2023ics0167947322002365.html)).
#' The `"continuous_nonparametric"` method does not assume anything about
#' the distributions of the comparison vectors. It directly
#' estimates the density ratio between the matches and the non-matches, using
#' the Kullback-Leibler Importance Estimation Procedure (KLIEP).
#' For details see [Sugiyama et al. (2008)](https://link.springer.com/article/10.1007/s10463-008-0197-x).
#'
#' @return
#' Returns a list containing:\cr
#' \itemize{
#' \item{`b_vars` -- a character vector of variables used for the `"binary"` method (with the prefix `"gamma_"`),}
#' \item{`cpar_vars` -- a character vector of variables used for the `"continuous_parametric"` method (with the prefix `"gamma_"`),}
#' \item{`cnonpar_vars` -- a character vector of variables used for the `"continuous_nonparametric"` method (with the prefix `"gamma_"`),}
#' \item{`b_params` -- parameters estimated using the `"binary"` method,}
#' \item{`cpar_params` -- parameters estimated using the `"continuous_parametric"` method,}
#' \item{`cnonpar_params` -- probability of exact matching estimated using the `"continuous_nonparametric"` method,}
#' \item{`ratio_kliep` -- a result of the \link[densityratio]{kliep} function,}
#' \item{`ratio_kliep_list` -- an object containing the results of the \link[densityratio]{kliep} function,}
#' \item{`ml_model` -- here `NULL`,}
#' \item{`pi_est` -- a prior probability of matching,}
#' \item{`match_prop` -- proportion of matches in the smaller dataset,}
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
#'   "name" = c("James", "Emma", "William", "Olivia", "Thomas",
#'   "Sophie", "Harry", "Amelia", "George", "Isabella"),
#'   "surname" = c("Smith", "Johnson", "Brown", "Taylor", "Wilson",
#'   "Davis", "Clark", "Harris", "Lewis", "Walker")
#' )
#'  df_2 <- data.frame(
#'   "name" = c("James", "Ema", "Wimliam", "Olivia", "Charlotte",
#'   "Henry", "Lucy", "Edward", "Alice", "Jack"),
#'   "surname" = c("Smith", "Johnson", "Bron", "Tailor", "Moore",
#'   "Evans", "Hall", "Wright", "Green", "King")
#' )
#' comparators <- list("name" = jarowinkler_complement(),
#'                     "surname" = jarowinkler_complement())
#' matches <- data.frame("a" = 1:4, "b" = 1:4)
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
    nonpar_hurdle = TRUE,
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
    prob_ratio <- "1"
  }

  # if (any(methods %in% c("binary", "continuous_parametric"))) {

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

  # }

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
  cnonpar_params <- NULL
  ratio_kliep <- NULL
  ratio_kliep_list <- NULL

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

      if (nonpar_hurdle) {

        p_0_M_cnonpar <- p_0_formula(M_cnonpar)
        p_0_U_cnonpar <- p_0_formula(Omega_cnonpar)

        ratio_kliep_list <- lapply(cnonpar_vars, function(x) {
          gamma_M <- M_cnonpar[[x]]
          gamma_M <- data.table::data.table(gamma_M[gamma_M > 0])
          names(gamma_M) <- x
          gamma_Omega <- Omega_cnonpar[[x]]
          gamma_Omega <- data.table::data.table(gamma_Omega[gamma_Omega > 0])
          names(gamma_Omega) <- x
          if (length(gamma_M) > 0) {
            ratio_plus <- do.call(
              densityratio::kliep,
              c(list(
                df_numerator = gamma_M,
                df_denominator = gamma_Omega
              ),
              controls_kliep)
            )
          } else {
            NULL
          }
        })
        names(ratio_kliep_list) <- cnonpar_vars

        cnonpar_params <- data.table(
          variable = cnonpar_vars,
          p_0_M_cnonpar = p_0_M_cnonpar,
          p_0_U_cnonpar = p_0_U_cnonpar
        )

      } else {

        ratio_kliep <- do.call(
          densityratio::kliep,
          c(list(
            df_numerator = M_cnonpar,
            df_denominator = Omega_cnonpar
          ),
          controls_kliep)
        )

      }

    } else if (prob_ratio == "2") {
      U_cnonpar <- U[, cnonpar_vars, with = FALSE]

      if (nonpar_hurdle) {

        p_0_M_cnonpar <- p_0_formula(M_cnonpar)
        p_0_U_cnonpar <- p_0_formula(U_cnonpar)

        ratio_kliep_list <- lapply(cnonpar_vars, function(x) {
          gamma_M <- M_cnonpar[[x]]
          gamma_M <- data.table::data.table(gamma_M[gamma_M > 0])
          names(gamma_M) <- x
          gamma_U <- U_cnonpar[[x]]
          gamma_U <- data.table::data.table(gamma_U[gamma_U > 0])
          names(gamma_U) <- x
          if (length(gamma_M) > 0) {
            ratio_plus <- do.call(
              densityratio::kliep,
              c(list(
                df_numerator = gamma_M,
                df_denominator = gamma_U
              ),
              controls_kliep)
            )
          } else {
            NULL
          }
        })
        names(ratio_kliep_list) <- cnonpar_vars

        cnonpar_params <- data.table(
          variable = cnonpar_vars,
          p_0_M_cnonpar = p_0_M_cnonpar,
          p_0_U_cnonpar = p_0_U_cnonpar
        )

      } else {

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

  }

  structure(
    list(
      b_vars = if (is.null(b_vars)) NULL else b_vars,
      cpar_vars = if (is.null(cpar_vars)) NULL else cpar_vars,
      cnonpar_vars = if (is.null(cnonpar_vars)) NULL else cnonpar_vars,
      b_params = if (is.null(b_params)) NULL else b_params,
      cpar_params = if (is.null(cpar_params)) NULL else cpar_params,
      cnonpar_params = if (is.null(cnonpar_params)) NULL else cnonpar_params,
      ratio_kliep = if (is.null(ratio_kliep)) NULL else ratio_kliep,
      ratio_kliep_list = if (is.null(ratio_kliep_list)) NULL else ratio_kliep_list,
      ml_model = NULL,
      pi_est = pi_est,
      match_prop = vectors$match_prop,
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
#' @details
#' The `custom_rec_lin_model` function creates a custom record linkage model,
#' based on known matches and non-matches (which might later serve as a classifier
#' for pairs outside training data). The procedure of creating a custom model
#' based on training data is as follows.
#' \enumerate{
#' \item{Use the `comparison_vectors` function to compare pairs of records.}
#' \item{Train a machine learning classifier using the `Omega` element
#' of the output of the `comparison_vectors` function. The classifier should
#' predict the probability of matching based on a given vector.}
#' \item{Use the `custom_rec_lin_model` function with
#' appropriate arguments.}
#' }
#'
#' @return
#' Returns a list containing:\cr
#' \itemize{
#' \item{`b_vars` -- here `NULL`,}
#' \item{`cpar_vars` -- here `NULL`,}
#' \item{`cnonpar_vars` -- here `NULL`,}
#' \item{`b_params` -- here `NULL`,}
#' \item{`cpar_params` -- here `NULL`,}
#' \item{`cnonpar_params` -- here `NULL`,}
#' \item{`ratio_kliep` -- here `NULL`,}
#' \item{`ratio_kliep_list` -- here `NULL`,}
#' \item{`ml_model` -- ML model used for creating the record linkage model,}
#' \item{`pi_est` -- a prior probability of matching,}
#' \item{`match_prop` -- proportion of matches in the smaller dataset,}
#' \item{`variables` -- a character vector of key variables used for comparison,}
#' \item{`comparators` -- a list of functions used to compare pairs of records,}
#' \item{`methods` -- here `NULL`.}
#' }
#'
#' @examples
#' if (requireNamespace("xgboost", quietly = TRUE)) {
#'   df_1 <- data.frame(
#'     "name" = c("James", "Emma", "William", "Olivia", "Thomas",
#'     "Sophie", "Harry", "Amelia", "George", "Isabella"),
#'     "surname" = c("Smith", "Johnson", "Brown", "Taylor", "Wilson",
#'     "Davis", "Clark", "Harris", "Lewis", "Walker")
#'   )
#'   df_2 <- data.frame(
#'     "name" = c("James", "Ema", "Wimliam", "Olivia", "Charlotte",
#'     "Henry", "Lucy", "Edward", "Alice", "Jack"),
#'     "surname" = c("Smith", "Johnson", "Bron", "Tailor", "Moore",
#'     "Evans", "Hall", "Wright", "Green", "King")
#'   )
#'   comparators <- list("name" = jarowinkler_complement(),
#'                       "surname" = jarowinkler_complement())
#'   matches <- data.frame("a" = 1:4, "b" = 1:4)
#'   vectors <- comparison_vectors(A = df_1, B = df_2, variables = c("name", "surname"),
#'                                comparators = comparators, matches = matches)
#'   train_data <- xgboost::xgb.DMatrix(
#'     data = as.matrix(vectors$Omega[, c("gamma_name", "gamma_surname")]),
#'     label = vectors$Omega$match
#'   )
#'   params <- list(objective = "binary:logistic",
#'                  eval_metric = "logloss")
#'   model_xgb <- xgboost::xgboost(data = train_data, params = params,
#'                                 nrounds = 100, verbose = 0)
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
      cnonpar_params = NULL,
      ratio_kliep = NULL,
      ratio_kliep_list = NULL,
      ml_model = ml_model,
      pi_est = pi_est,
      match_prop = vectors$match_prop,
      variables = variables,
      comparators = vectors$comparators,
      methods = NULL
    ),
    class = "rec_lin_model"
  )

}

