#' @import data.table
#' @importFrom stats dbinom
#' @importFrom FixedPoint FixedPoint
#'
#' @title Predict matches based on a given record linkage model
#'
#' @export
predict.rec_lin_model <- function(object,
                                  newdata_A,
                                  newdata_B,
                                  set_construction = c("size", "flr"),
                                  fixed_method = "Newton",
                                  target_flr = 0.05,
                                  tol = 10^(-3),
                                  max_iter = 1000,
                                  ...) {

  stopifnot("`newdata_A` is required for predictions." =
              !missing(newdata_A))

  stopifnot("`newdata_B` is required for predictions." =
              !missing(newdata_B))

  stopifnot("`set_construction` should be `size` or `flr`." =
              set_construction %in% c("size", "flr"))

  if (missing(set_construction)) set_construction <- "size"

  vectors <- comparison_vectors(A = newdata_A,
                                B = newdata_B,
                                variables = object$variables,
                                comparators = object$comparators)
  Omega <- vectors$Omega

  n <- NROW(Omega)

  if (!is.null(object$ml_model)) {

    # predicted_probs <- do.call(
    #   predict,
    #   c(list(
    #     object$ml_model,
    #     Omega
    #   ),
    #   controls_ml_predictions)
    # )
    ml_model <- object$ml_model
    predicted_probs <- predict(ml_model,
                               Omega,
                               ...)
    Omega[, "ratio" := predicted_probs]
    # TODO

  } else {

    Omega[, "ratio" := 1]

    if (!is.null(object$binary_variables)) {

      binary_variables <- object$binary_variables
      binary_params <- object$binary_params
      Omega_binary <- Omega[, ..binary_variables]
      binary_numerator_list <- lapply(binary_variables,
                                      function(col) {
                                        stats::dbinom(x = Omega_binary[[col]],
                                                      size = 1,
                                                      prob = as.numeric(binary_params[variable == col, "theta"]))
                                      })
      binary_numerator <- Reduce(`*`, binary_numerator_list)
      binary_denominator_list <- lapply(binary_variables,
                                        function(col) {
                                          stats::dbinom(x = Omega_binary[[col]],
                                                        size = 1,
                                                        prob = as.numeric(binary_params[variable == col, "eta"]))
                                        })
      binary_denominator <- Reduce(`*`, binary_denominator_list)
      Omega[, ratio := ratio * binary_numerator / binary_denominator]

    }

    if (!is.null(object$continuous_parametric_variables)) {

      continuous_parametric_variables <- object$continuous_parametric_variables
      continuous_parametric_params <- object$continuous_parametric_params
      Omega_continuous_parametric <- Omega[, ..continuous_parametric_variables]
      continuous_parametric_numerator_list <- lapply(continuous_parametric_variables,
                                                     function(col) {
                                                       hurdle_gamma_density(x = Omega_continuous_parametric[[col]],
                                                                            p_0 = as.numeric(continuous_parametric_params[variable == col, "p_0_M"]),
                                                                            alpha = as.numeric(continuous_parametric_params[variable == col, "alpha_M"]),
                                                                            beta = as.numeric(continuous_parametric_params[variable == col, "beta_M"]))
                                                     })
      continuous_parametric_numerator <- Reduce(`*`, continuous_parametric_numerator_list)
      continuous_parametric_denominator_list <- lapply(continuous_parametric_variables,
                                                       function(col) {
                                                         hurdle_gamma_density(x = Omega_continuous_parametric[[col]],
                                                                              p_0 = as.numeric(continuous_parametric_params[variable == col, "p_0_U"]),
                                                                              alpha = as.numeric(continuous_parametric_params[variable == col, "alpha_U"]),
                                                                              beta = as.numeric(continuous_parametric_params[variable == col, "beta_U"]))
                                                       })
      continuous_parametric_denominator <- Reduce(`*`, continuous_parametric_denominator_list)
      Omega[, ratio := ratio * continuous_parametric_numerator / continuous_parametric_denominator]

    }

    if (!is.null(object$continuous_nonparametric_variables)) {

      continuous_nonparametric_variables <- object$continuous_nonparametric_variables
      ratio_kliep <- object$ratio_kliep
      Omega_continuous_nonparametric <- Omega[, ..continuous_nonparametric_variables]
      Omega[, ratio := ratio * predict(ratio_kliep, Omega_continuous_nonparametric)]

    }

  }

  n_M_start <- NROW(merge(newdata_A, newdata_B, by = object$variables, all = FALSE))
  fun_n_M <- fixed_n_M(n = n, ratio_gamma = Omega$ratio)
  n_M_est <- min(FixedPoint::FixedPoint(Function = fun_n_M, # to think
                                        Inputs = n_M_start,
                                        Method = fixed_method)$FixedPoint,
                 min(NROW(newdata_A), NROW(newdata_B)))
  n_M_est <- max(n_M_est, 0)

  if (set_construction == "size") {

    # n_M_start <- NROW(merge(newdata_A, newdata_B, by = object$variables, all = FALSE))
    # fun_n_M <- fixed_n_M(n = n, ratio_gamma = Omega$ratio)
    # n_M_est <- min(FixedPoint::FixedPoint(Function = fun_n_M, # to think
    #                                       Inputs = n_M_start,
    #                                       Method = fixed_method)$FixedPoint,
    #                min(NROW(newdata_A), NROW(newdata_B)))
    # n_M_est <- max(n_M_est, 0)

    Omega <- Omega[order(-ratio), ]
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

    g_est <- pmin(NROW(M_est) * M_est$ratio / (NROW(M_est) * (M_est$ratio - 1) + n), 1)
    flr_est <- 1 / NROW(M_est) * sum(1 - g_est)

    it <- NULL

  } else if (set_construction == "flr") {

    Omega <- Omega[order(-ratio), ]

    min_treshold <- min(Omega$ratio)
    max_treshold <- max(Omega$ratio)
    treshold <- (min_treshold + max_treshold) / 2

    it <- 0

    while (it < max_iter) {

      M_est <- Omega[ratio >= treshold, ]
      g_est <- pmin(NROW(M_est) * M_est$ratio / (NROW(M_est) * (M_est$ratio - 1) + n), 1)
      flr_est <- 1 / NROW(M_est) * sum(1 - g_est)

      if (abs(flr_est - target_flr) <= tol) {

        break

      } else if (flr_est < target_flr) {

        max_treshold <- treshold
        treshold <- (min_treshold + max_treshold) / 2

      } else {

        min_treshold <- treshold
        treshold <- (min_treshold + max_treshold) / 2

      }

      it <- it + 1

    }

  }

  mmr_est <- 1 - sum(g_est / n_M_est)
  M_est <- M_est[, c("a", "b")]

  structure(
    list(
      M_est = M_est,
      set_construction = set_construction,
      flr_est = flr_est,
      mmr_est = mmr_est,
      it = if (is.null(it)) NULL else it
    ),
    class = "rec_lin_predictions"
  )

}
