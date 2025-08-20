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
#' Predicts matches between records in two datasets based on a given record linkage model.
#'
#' @param object A `rec_lin_model` object from the `train_rec_lin` or `custom_rec_lin_model` functions.
#' @param newdata_A A duplicate-free `data.frame` or `data.table`.
#' @param newdata_B A duplicate-free `data.frame` or `data.table`.
#' @param set_construction A method for constructing the predicted set of matches (`"size"` or `"flr"`).
#' @param fixed_method A method for solving fixed-point equations using the \link[FixedPoint]{FixedPoint} function
#' (used only if `set_construction == "size"`).
#' @param target_flr A target false link rate (FLR) (used only if `set_construction == "flr"`).
#' @param tol Error tolerance in the bisection procedure (used only if `set_construction == "flr"`).
#' @param max_iter A maximum number of iterations for the bisection procedure (used only if `set_construction == "flr"`).
#' @param data_type Data type for predictions with a custom ML model (`"data.frame"`, `"data.table"` or `"matrix"`;
#' used only if `object` is from the `custom_rec_lin_model` function).
#' @param true_matches A `data.frame` or `data.table` indicating true matches.
#' @param ... Additional controls passed to the `predict` function for custom ML model
#' (used only if the `object` is from the `custom_rec_lin_model` function).
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
#' \item{`eval_metrics` -- metrics for quality assessment, if `true_matches` is provided,}
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
#'   "name" = c("John", "Emily", "Mark", "Anna", "David"),
#'   "surname" = c("Smith", "Johnson", "Taylor", "Williams", "Brown")
#' )
#' df_2 <- data.frame(
#'   "name" = c("John", "Emely", "Marc", "Michael"),
#'   "surname" = c("Smith", "Jonson", "Tailor", "Henderson")
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
#'
#' df_new_1 <- data.frame(
#'   "name" = c("Jame", "Lia", "Tomas", "Matthew", "Andrew"),
#'   "surname" = c("Wilsen", "Thomsson", "Davis", "Robinson", "Scott")
#' )
#' df_new_2 <- data.frame(
#'   "name" = c("James", "Leah", "Thomas", "Sophie", "Mathew", "Andrew"),
#'   "surname" = c("Wilson", "Thompson", "Davies", "Clarks", "Robins", "Scots")
#' )
#' predict(model, df_new_1, df_new_2, set_construction = "size")
#' @export
predict.rec_lin_model <- function(object,
                                  newdata_A,
                                  newdata_B,
                                  set_construction = c("size", "flr"),
                                  fixed_method = "Newton",
                                  target_flr = 0.05,
                                  tol = 10^(-3),
                                  max_iter = 50,
                                  data_type = c("data.frame", "data.table", "matrix"),
                                  true_matches = NULL,
                                  ...) {

  stopifnot("`newdata_A` is required for predictions." =
              !missing(newdata_A))

  stopifnot("`newdata_B` is required for predictions." =
              !missing(newdata_B))

  stopifnot("`set_construction` should be `size` or `flr`." =
              set_construction %in% c("size", "flr"))

  if (missing(set_construction)) set_construction <- "size"

  if (missing(data_type)) data_type <- "data.frame"

  vectors <- comparison_vectors(A = newdata_A,
                                B = newdata_B,
                                variables = object$variables,
                                comparators = object$comparators)
  Omega <- vectors$Omega

  n <- NROW(Omega)

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

    predicted_ratio <- predicted_probs * (1 - object$pi_est) / ((1 - predicted_probs) * object$pi_est)
    data.table::set(Omega, j = "ratio", value = predicted_ratio)

  } else {

    data.table::set(Omega, j = "ratio", value = 1)


    if (!is.null(object$b_vars)) {

      b_vars <- object$b_vars
      b_params <- object$b_params
      Omega_b <- Omega[, b_vars, with = FALSE]
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

    }

    if (!is.null(object$cpar_vars)) {

      cpar_vars <- object$cpar_vars
      cpar_params <- object$cpar_params
      if ("p_0_Omega" %in% colnames(cpar_params)) {
        data.table::setnames(cpar_params,
                             old = c("p_0_Omega", "alpha_Omega", "beta_Omega"),
                             new = c("p_0_U", "alpha_U", "beta_U"))
      }
      Omega_cpar <- Omega[, cpar_vars, with = FALSE]
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

    }

    if (!is.null(object$cnonpar_vars)) {

      cnonpar_vars <- object$cnonpar_vars
      ratio_kliep <- object$ratio_kliep
      Omega_cnonpar <- Omega[, cnonpar_vars, with = FALSE]
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * stats::predict(ratio_kliep, Omega_cnonpar))

    }

  }

  # to think
  n_M_start <- min(NROW(newdata_A), NROW(newdata_B))
  fun_n_M <- fixed_n_M(n = n, ratio_gamma = Omega[["ratio"]])
  n_M_original <- FixedPoint::FixedPoint(Function = fun_n_M,
                                         Inputs = n_M_start,
                                         Method = fixed_method)$FixedPoint
  # n_M_est <- min(FixedPoint::FixedPoint(Function = fun_n_M,
  #                                       Inputs = n_M_start,
  #                                       Method = fixed_method)$FixedPoint,
  #                min(NROW(newdata_A), NROW(newdata_B)))
  n_M_est <- min(n_M_original, n_M_start)

  n_M_est <- max(n_M_est, 0)
  n_M_est <- round(n_M_est)

  if (set_construction == "size") {

    Omega <- Omega[order(-get("ratio")), ]
    M_est <- data.table("a" = numeric(), "b" = numeric(), "ratio" = numeric())
    used_a <- c()
    used_b <- c()

    for (i in 1:NROW(Omega)) {

      current_a <- Omega$a[i]
      current_b <- Omega$b[i]
      if (!(current_a %in% used_a) && !(current_b %in% used_b)) {
        M_est <- rbind(M_est, Omega[i, c("a", "b", "ratio")])
        used_a <- c(used_a, current_a)
        used_b <- c(used_b, current_b)
      }

    }

    M_est <- head(M_est, round(n_M_est))
    g_est <- pmin(NROW(M_est) * M_est$ratio / (NROW(M_est) * (M_est$ratio - 1) + n), 1)
    flr_est <- 1 / NROW(M_est) * sum(1 - g_est)

    iter <- NULL

    mmr_est <- 1 - sum(g_est / n_M_est)

  } else if (set_construction == "flr") {

    Omega <- Omega[order(-get("ratio")), ]

    min_treshold <- min(Omega$ratio)
    max_treshold <- max(Omega$ratio)
    treshold <- (min_treshold + max_treshold) / 2

    iter <- 0

    while (iter < max_iter) {

      M_est <- Omega[get("ratio") >= treshold, ]
      M_est$g_est <- pmin(NROW(M_est) * M_est$ratio / (NROW(M_est) * (M_est$ratio - 1) + n), 1)
      flr_est <- 1 / NROW(M_est) * sum(1 - M_est$g_est)

      if (abs(flr_est - target_flr) <= tol) {

        break

      } else if (flr_est < target_flr) {

        max_treshold <- treshold
        treshold <- (min_treshold + max_treshold) / 2

      } else {

        min_treshold <- treshold
        treshold <- (min_treshold + max_treshold) / 2

      }

      iter <- iter + 1

    }

    g_est_to_mmr <- pmin(n_M_original * M_est$ratio / (n_M_original * (M_est$ratio - 1) + n), 1)

    # mmr_est <- 1 - sum(M_est$g_est / n_M_original)
    mmr_est <- 1 - sum(g_est_to_mmr / n_M_original) # to think

  }

  M_est <- M_est[, c("a", "b", "ratio")]

  if (!is.null(true_matches)) {

    data.table::setDT(true_matches)

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
