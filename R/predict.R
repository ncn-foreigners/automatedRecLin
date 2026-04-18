#' @import data.table
#' @importFrom stats dbinom
#' @importFrom stats predict
#' @importFrom FixedPoint FixedPoint
#' @importFrom utils head
#'
#' @title Predict Matches Based on a Given Record Linkage Model
#'
#' @author Adam Struzik
#'
#' @description
#' Predicts matches between records in two datasets based on a given record linkage model,
#' using the maximum entropy classification (MEC) algorithm
#' (see [Lee et al. (2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm)).
#'
#' @param object A `rec_lin_model` object from the `train_rec_lin` or `custom_rec_lin_model` functions.
#' @param newdata_A A duplicate-free `data.frame` or `data.table`.
#' @param newdata_B A duplicate-free `data.frame` or `data.table`.
#' @param duplicates_in_A Logical indicating whether to allow `A` to have duplicate records.
#' @param set_construction A method for constructing the predicted set of matches (`"size"`, `"flr"` or `"mmr"`).
#' @param fixed_method A method for solving fixed-point equations using the \link[FixedPoint]{FixedPoint} function.
#' @param target_rate A target false link rate (FLR) or missing match rate (MMR)
#' (used only if `set_construction == "flr"` or `set_construction == "mmr"`).
#' @param tol Error tolerance in the bisection procedure
#' (used only if `set_construction == "flr"` or `set_construction == "mmr"`).
#' @param max_iter A maximum number of iterations for the bisection procedure
#' (used only if `set_construction == "flr"` or `set_construction == "mmr"`).
#' @param data_type Data type for predictions with a custom ML model (`"data.frame"`, `"data.table"` or `"matrix"`;
#' used only if `object` is from the `custom_rec_lin_model` function).
#' @param true_matches A `data.frame` or `data.table` indicating true matches.
#' @param ... Additional controls passed to the `predict` function for custom ML model
#' (used only if the `object` is from the `custom_rec_lin_model` function).
#'
#' @details
#' The `predict` function estimates the probability/density ratio
#' between matches and non-matches for pairs in given
#' datasets, based on a model obtained using the
#' `train_rec_lin` or `custom_rec_lin_model` functions.
#' Then, it estimates the number of matches and
#' returns the predicted matches, using the maximum
#' entropy classification (MEC) algorithm
#' (see [Lee et al. (2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm)).
#'
#' The `predict` function allows the construction of the predicted set
#' of matches using its estimated size or the bisection procedure,
#' described in [Lee et al. (2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm),
#' based on a target False Link Rate (FLR)
#' or missing match rate (MMR). To use the second option, set `set_construction = "flr"`
#' or `set_construction = "mmr"` and
#' specify a target error rate using the `target_rate` argument.
#'
#' By default, the function assumes that the datasets `newdata_A` and `newdata_B`
#' contain no duplicate records. This assumption
#' might be relaxed by allowing `newdata_A` to have duplicates. To do so,
#' set `duplicates_in_A = TRUE`.
#'
#' @return
#' Returns a list containing:\cr
#' \itemize{
#' \item{`M_est` -- a `data.table` with predicted matches,}
#' \item{`set_construction` -- a method for constructing the predicted set of matches,}
#' \item{`n_M_est` -- estimated classification set size,}
#' \item{`flr_est` -- estimated false link rate (FLR),}
#' \item{`mmr_est` -- estimated missing match rate (MMR),}
#' \item{`iter` -- the number of iterations in the bisection procedure,}
#' \item{`eval_metrics` -- standard metrics for quality assessment, if `true_matches` is provided,}
#' \item{`confusion` -- confusion matrix, if `true_matches` is provided.}
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
#'
#' df_new_1 <- data.frame(
#'   "name" = c("John", "Emily", "Mark", "Anna", "David"),
#'   "surname" = c("Smith", "Johnson", "Taylor", "Williams", "Brown")
#' )
#' df_new_2 <- data.frame(
#'   "name" = c("John", "Emely", "Mark", "Michael"),
#'   "surname" = c("Smitth", "Johnson", "Tailor", "Henders")
#' )
#' predict(model, df_new_1, df_new_2)
#' @export
predict.rec_lin_model <- function(object,
                                  newdata_A,
                                  newdata_B,
                                  duplicates_in_A = FALSE,
                                  set_construction = c("size", "flr", "mmr"),
                                  fixed_method = "Newton",
                                  target_rate = 0.03,
                                  tol = 0.005,
                                  max_iter = 50,
                                  data_type = c("data.frame", "data.table", "matrix"),
                                  true_matches = NULL,
                                  ...) {

  stopifnot("`newdata_A` is required for predictions." =
              !missing(newdata_A))

  stopifnot("`newdata_B` is required for predictions." =
              !missing(newdata_B))

  if (missing(set_construction)) set_construction <- "size"
  if (missing(data_type)) data_type <- "data.frame"
  validate_choice(set_construction, c("size", "flr", "mmr"), "set_construction")
  validate_choice(data_type, c("data.frame", "data.table", "matrix"), "data_type")
  true_matches <- sanitize_true_matches(
    true_matches = true_matches,
    n_A = nrow(newdata_A),
    n_B = nrow(newdata_B),
    arg_name = "true_matches"
  )

  vectors <- comparison_vectors(A = newdata_A,
                                B = newdata_B,
                                variables = object$variables,
                                comparators = object$comparators)
  Omega <- vectors$Omega

  n <- NROW(Omega)
  prob_est <- object$match_prop / max(NROW(newdata_A), NROW(newdata_B))

  if (!is.null(object$ml_model)) {

    ml_model <- object$ml_model
    gamma_variables <- paste0("gamma_", vectors$variables)

    if (data_type == "data.frame") {

      Omega_temp <- as.data.frame(Omega[, gamma_variables, with = FALSE])

    } else if (data_type == "matrix") {

      Omega_temp <- as.matrix(Omega[, gamma_variables, with = FALSE])

    } else {

      Omega_temp <- Omega

    }

    controls_ml_predictions <- list(...)

    predicted_probs <- do.call(
      predict,
      c(list(
        ml_model,
        Omega_temp
      ),
      controls_ml_predictions)
    )

    predicted_ratio <- predicted_probs * (1 - prob_est) / ((1 - predicted_probs) * prob_est)
    data.table::set(Omega, j = "ratio", value = predicted_ratio)

  } else {

    data.table::set(Omega, j = "ratio", value = 1)


    if (length(object$b_vars) > 0) {

      b_vars <- object$b_vars
      b_params <- align_parameter_table(object$b_params, b_vars)
      Omega_b <- Omega[, b_vars, with = FALSE]
      data.table::set(
        Omega,
        j = "ratio",
        value = Omega[["ratio"]] * bernoulli_ratio(Omega_b, b_params$theta, b_params$eta)
      )

    }

    if (length(object$cpar_vars) > 0) {

      cpar_vars <- object$cpar_vars
      cpar_params <- data.table::copy(align_parameter_table(object$cpar_params, cpar_vars))
      if (object$prob_ratio == "1") {
        p_0_denominator <- cpar_params$p_0_Omega
        alpha_denominator <- cpar_params$alpha_Omega
        beta_denominator <- cpar_params$beta_Omega
      } else {
        p_0_denominator <- cpar_params$p_0_U
        alpha_denominator <- cpar_params$alpha_U
        beta_denominator <- cpar_params$beta_U
      }
      Omega_cpar <- Omega[, cpar_vars, with = FALSE]
      data.table::set(
        Omega,
        j = "ratio",
        value = Omega[["ratio"]] * hurdle_gamma_ratio(
          Omega_cpar,
          cpar_params$p_0_M,
          cpar_params$alpha_M,
          cpar_params$beta_M,
          p_0_denominator,
          alpha_denominator,
          beta_denominator
        )
      )

    }

    if (length(object$cnonpar_vars) > 0) {

      cnonpar_vars <- object$cnonpar_vars
      Omega_cnonpar <- Omega[, cnonpar_vars, with = FALSE]

      if (!is.null(object$ratio_kliep_list)) {

        ratio_kliep_list <- object$ratio_kliep_list
        missing_models <- missing_kliep_models(ratio_kliep_list)
        if (length(missing_models) > 0L) {
          warn_kliep_issue("predict()", missing_models, "using only the hurdle mass term for those variables.")
        }
        p_0_M_cnonpar <- object$cnonpar_params[["p_0_M_cnonpar"]]
        p_0_U_cnonpar <- object$cnonpar_params[["p_0_U_cnonpar"]]
        names(p_0_M_cnonpar) <- cnonpar_vars
        names(p_0_U_cnonpar) <- cnonpar_vars

        ratio_kliep <- kliep_hurdle_ratio(
          Omega_cnonpar,
          cnonpar_vars,
          p_0_M_cnonpar,
          p_0_U_cnonpar,
          ratio_kliep_list
        )
        data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * ratio_kliep)

      } else {

        ratio_kliep <- object$ratio_kliep
        data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * stats::predict(ratio_kliep, Omega_cnonpar))

      }

    }

  }

  n_M_start <- min(NROW(newdata_A), NROW(newdata_B))
  if (object$prob_ratio == "1" && is.null(object$ml_model)) {
    # Under ratio-I, r_q = m / q and the posterior match probability is pi * r_q.
    g_est <- pmin(prob_est * Omega[["ratio"]], 1)
    n_M_est <- min(sum(g_est), n_M_start)
  } else {
    fun_n_M <- fixed_n_M(n = n, ratio_gamma = Omega[["ratio"]])
    n_M_original <- FixedPoint::FixedPoint(Function = fun_n_M,
                                           Inputs = n_M_start,
                                           Method = fixed_method)$FixedPoint
    n_M_est <- min(n_M_original, n_M_start)
    n_M_est <- max(n_M_est, 0)
    n_M_est <- round(n_M_est)

    g_est <- pmin(n_M_est * Omega$ratio / (n_M_est * (Omega$ratio - 1) + n), 1)
  }
  selection_summary <- summarize_mec_selection(
    a = Omega[["a"]],
    b = Omega[["b"]],
    ratio = Omega[["ratio"]],
    g_est = g_est,
    n_M_est = n_M_est,
    duplicates_in_A = duplicates_in_A,
    set_construction = set_construction,
    target_rate = target_rate,
    tol = tol,
    max_iter = max_iter
  )
  M_est <- Omega[selection_summary$selected_idx, c("a", "b", "ratio")]
  flr_est <- selection_summary$flr_est
  mmr_est <- selection_summary$mmr_est
  iter <- selection_summary$iter

  if (!is.null(true_matches)) {
    eval <- evaluation(M_est, true_matches, n)
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
      M_est = M_est,
      set_construction = set_construction,
      n_M_est = n_M_est,
      flr_est = flr_est,
      mmr_est = mmr_est,
      iter = if (is.null(iter)) NULL else iter,
      eval_metrics = if (is.null(true_matches)) NULL else eval_metrics,
      confusion = if (is.null(true_matches)) NULL else confusion
    ),
    class = "rec_lin_predictions"
  )

}
