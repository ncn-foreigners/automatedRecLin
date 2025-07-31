#' @method print comparison_vectors
#' @exportS3Method
print.comparison_vectors <- function(x, ...) {

  cat("Comparison based on the following variables: ", paste(x$variables, collapse = ", "), ".\n", sep = "")
  if (!("match" %in% colnames(x$Omega))) {
    cat("True matches are unknown.\n")
  }
  cat("========================================================\n")
  print(x$Omega)

}

#' @method print rec_lin_model
#' @exportS3Method
print.rec_lin_model <- function(x, ...) {

  cat("Record linkage model based on the following variables: ", paste(x$variables, collapse = ", "), ".\n", sep = "")

  if (!is.null(x$ml_model)) {
    cat("A custom ML model was used.\n")
  }

  cat("The prior probability of matching is ", x$pi_est, ".\n", sep = "")

  if (!is.null(x$binary_variables)) {
    cat("========================================================\n")
    cat("Variables selected for the binary method: ", paste(substring(x$binary_variables, 7), collapse = ", "), ".\n", sep = "")
    cat("Estimated parameters for the binary method:\n")
    print(x$binary_params)
  }

  if (!is.null(x$continuous_parametric_variables)) {
    cat("========================================================\n")
    cat("Variables selected for the continuous parametric method: ", paste(substring(x$continuous_parametric_variables, 7), collapse = ", "), ".\n", sep = "")
    cat("Estimated parameters for the continuous parametric method:\n")
    print(x$binary_params)
  }

  if (!is.null(x$continuous_nonparametric_variables)) {
    cat("========================================================\n")
    cat("Variables selected for the continuous nonparametric method: ", paste(substring(x$continuous_nonparametric_variables, 7), collapse = ", "), ".\n", sep = "")
  }

}

#' @method print rec_lin_predictions
#' @exportS3Method
print.rec_lin_predictions <- function(x, ...) {

  if (NROW(x$M_est) == 0) {
    cat("No matches were predicted.\n")
  } else {
    cat("Predicted matches:\n")
    print(x$M_est)
    cat("========================================================\n")
    if (x$set_construction == "size") {
      cat("The construction of the classification set was based on estimates of its size.\n")
    } else if (x$set_construction == "flr") {
      cat("The construction of the classification set was based on the target false link rate (FLR).\n")
    }
    cat("Estimated false link rate (FLR): ", x$flr_est, ".\n", sep = "")
    cat("Estimated missing match rate (MMR): ", x$mmr_est, ".\n", sep = "")
  }

}
