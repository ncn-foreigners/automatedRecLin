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
                                  controls_ml_predictions = list()) {

  stopifnot("`newdata_A` is required for predictions." =
              !missing(newdata_A))

  stopifnot("`newdata_B` is required for predictions." =
              !missing(newdata_B))

  if (missing(set_construction)) set_construction <- "size"

  vectors <- comparison_vectors(A = newdata_A,
                                B = newdata_B,
                                variables = object$variables,
                                comparators = object$comparators)
  Omega <- vectors$Omega

  n <- NROW(Omega)

  if (!is.null(object$ml_model)) {

    predicted_probs <- do.call(
      predict,
      c(list(
        object$ml_model,
        Omega
      ),
      controls_ml_predictions)
    )
    Omega[, "ratio" := predicted_probs]

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

  if (set_construction == "size") {

    n_M_start <- NROW(merge(newdata_A, newdata_B, by = object$variables, all = FALSE))
    fun_n_M <- fixed_n_M(n = n, ratio_gamma = Omega$ratio)
    n_M_est <- min(round(FixedPoint::FixedPoint(Function = fun_n_M, # to think
                                                Inputs = n_M_start,
                                                Method = fixed_method)$FixedPoint),
                   min(NROW(newdata_A), NROW(newdata_B)))

    Omega <- Omega[order(-ratio), ]
    M_est <- data.table("a" = numeric(), "b" = numeric())
    used_a <- c()
    used_b <- c()

    for (i in 1:NROW(Omega)) {

      current_a <- Omega$a[i]
      current_b <- Omega$b[i]
      if (!(current_a %in% used_a) && !(current_b %in% used_b)) {
        M_est <- rbind(M_est, Omega[i, c("a", "b")])
        used_a <- c(used_a, current_a)
        used_b <- c(used_b, current_b)
      }

    }

    return(head(M_est, n_M_est))

  }

}
