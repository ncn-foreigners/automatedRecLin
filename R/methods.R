#' @importFrom utils head
#' @method print comparison_vectors
#' @exportS3Method
print.comparison_vectors <- function(x, ...) {

  cat("Comparison based on the following variables: ", paste(x$variables, collapse = ", "), ".\n", sep = "")
  if (!("match" %in% colnames(x$Omega))) {
    cat("True matches are unknown.\n")
  }
  cat("========================================================\n")
  print(head(x$Omega))

}

#' @importFrom utils head
#' @method print rec_lin_model
#' @exportS3Method
print.rec_lin_model <- function(x, ...) {

  cat("Record linkage model based on the following variables: ", paste(x$variables, collapse = ", "), ".\n", sep = "")

  if (!is.null(x$ml_model)) {
    cat("A custom ML model was used.\n")
  }

  cat("The prior probability of matching is ", x$pi_est, ".\n", sep = "")

  if (!is.null(x$b_vars)) {
    cat("========================================================\n")
    cat("Variables selected for the binary method: ", paste(substring(x$b_vars, 7), collapse = ", "), ".\n", sep = "")
    cat("Estimated parameters for the binary method:\n")
    print(x$b_params)
  }

  if (!is.null(x$cpar_vars)) {
    cat("========================================================\n")
    cat("Variables selected for the continuous parametric method: ", paste(substring(x$cpar_vars, 7), collapse = ", "), ".\n", sep = "")
    cat("Estimated parameters for the continuous parametric method:\n")
    print(x$cpar_params)
  }

  if (!is.null(x$cnonpar_vars)) {
    cat("========================================================\n")
    cat("Variables selected for the continuous nonparametric method: ", paste(substring(x$cnonpar_vars, 7), collapse = ", "), ".\n", sep = "")
  }

  cat("========================================================\n")
  cat("Probability/density ratio type: ", x$prob_ratio, ".\n", sep = "")

}

#' @import data.table
#' @importFrom utils head
#' @method print rec_lin_predictions
#' @exportS3Method
print.rec_lin_predictions <- function(x, ...) {

  if (NROW(x$M_est) == 0) {
    cat("No matches were predicted.\n")
  } else {
    # cat("Predicted matches:\n")
    # print(x$M_est)
    # cat("========================================================\n")
    M_est <- data.table::copy(x$M_est)
    data.table::set(M_est, j = "ratio / 1000", value = M_est[["ratio"]] / 1000)
    data.table::set(M_est, j = "ratio", value = NULL)
    M_est_head <- head(M_est, 6)
    cat("The algorithm predicted", NROW(M_est), "matches.\n")
    cat("The first", NROW(M_est_head), "predicted matches are:\n")
    print(M_est_head)
    cat("========================================================\n")
    if (x$set_construction == "size") {
      cat("The construction of the classification set was based on estimates of its size.\n")
    } else if (x$set_construction == "flr") {
      cat("The construction of the classification set was based on the target false link rate (FLR).\n")
      cat("The bisection procedure ended after", x$iter, "iterations.\n")
    }
    cat("Estimated false link rate (FLR): ", sprintf("%.4f", x$flr_est * 100), " %.\n", sep = "")
    if (x$n_M_est != 0) {
      cat("Estimated missing match rate (MMR): ", sprintf("%.4f", x$mmr_est * 100), " %.\n", sep = "")
    } else {
      # cat("Estimated classification set size is 0.\n")
      cat("Missing match rate (MMR) cannot be estimated because the estimated classification set size is equal to 0.")
    }
  }

  if (!is.null(x$eval_metrics)) {
    cat("========================================================\n")
    cat("Evaluation metrics (presented in percentages):\n")
    eval_metrics <- as.numeric(sprintf("%.4f", x$eval_metrics * 100))
    names(eval_metrics) <- names(x$eval_metrics)
    print(eval_metrics)
    cat("Note that precision = 1 - flr, and fnr = mmr.")
  }

}

#' @import data.table
#' @importFrom utils head
#' @method print mec_rec_lin
#' @exportS3Method
print.mec_rec_lin <- function(x, ...) {

  cat("Record linkage based on the following variables: ", paste(x$variables, collapse = ", "), ".\n", sep = "")
  cat("========================================================\n")
  if (NROW(x$M_est) == 0) {
    cat("No matches were predicted.\n")
  } else {
    M_est <- data.table::copy(x$M_est)
    data.table::set(M_est, j = "ratio / 1000", value = M_est[["ratio"]] / 1000)
    data.table::set(M_est, j = "ratio", value = NULL)
    M_est_head <- head(M_est, 6)
    cat("The algorithm predicted", NROW(M_est), "matches.\n")
    cat("The first", NROW(M_est_head), "predicted matches are:\n")
    print(M_est_head)
    cat("========================================================\n")
  }
  if (x$set_construction == "size") {
    cat("The construction of the classification set was based on estimates of its size.\n")
  } else if (x$set_construction == "flr") {
  cat("The construction of the classification set was based on the target false link rate (FLR).\n")
  cat("The bisection procedure ended after", x$iter_bisection, "iterations.\n")
  }
  # cat("Estimated false link rate (FLR): ", x$flr_est, ".\n", sep = "")
  cat("Estimated false link rate (FLR): ", sprintf("%.4f", x$flr_est * 100), " %.\n", sep = "")
  if (!is.null(x$mmr_est)) {
    # cat("Estimated missing match rate (MMR): ", x$mmr_est, ".\n", sep = "")
    cat("Estimated missing match rate (MMR): ", sprintf("%.4f", x$mmr_est * 100), " %.\n", sep = "")
  }

  if (!is.null(x$b_vars)) {
    cat("========================================================\n")
    cat("Variables selected for the binary method: ", paste(substring(x$b_vars, 7), collapse = ", "), ".\n", sep = "")
    cat("Estimated parameters for the binary method:\n")
    print(x$b_params)
  }

  if (!is.null(x$cpar_vars)) {
    cat("========================================================\n")
    cat("Variables selected for the continuous parametric method: ", paste(substring(x$cpar_vars, 7), collapse = ", "), ".\n", sep = "")
    cat("Estimated parameters for the continuous parametric method:\n")
    print(x$cpar_params)
  }

  if (!is.null(x$cnonpar_vars)) {
    cat("========================================================\n")
    cat("Variables selected for the continuous nonparametric method: ", paste(substring(x$cnonpar_vars, 7), collapse = ", "), ".\n", sep = "")
  }

  if (!is.null(x$hm_params)) {
    cat("========================================================\n")
    cat("Variables selected for the hit-miss method: ", paste(substring(x$hm_vars, 7), collapse = ", "), ".\n", sep = "")
    cat("Estimated parameters for the hit-miss method:\n")
    print(x$hm_params)
  }

  if (!is.null(x$eval_metrics)) {
    cat("========================================================\n")
    cat("Evaluation metrics (presented in percentages):\n")
    eval_metrics <- as.numeric(sprintf("%.4f", x$eval_metrics * 100))
    names(eval_metrics) <- names(x$eval_metrics)
    print(eval_metrics)
    cat("Note that precision = 1 - flr, and fnr = mmr.")
  }

}
