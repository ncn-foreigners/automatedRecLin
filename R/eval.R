#' @title Evaluation for Record Linkage
#'
#' @noRd
evaluation <- function(pred_matches, true_matches, n) {

  pred_matches <- paste0(pred_matches[["a"]], "_", pred_matches[["b"]])
  true_matches <- paste0(true_matches[["a"]], "_", true_matches[["b"]])

  TP <- sum(pred_matches %in% true_matches)
  FP <- length(pred_matches) - TP
  TN <- n - length(true_matches) - FP
  FN <- sum(!(true_matches %in% pred_matches))

  list(
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN
  )
}

#' @title Evaluation Metrics for Record Linkage
#'
#' @noRd
get_metrics <- function(TP, FP, FN, TN) {

  precision <- if (TP + FP != 0) TP / (TP + FP) else 0
  fnr <- if (FN + TP != 0) FN / (FN + TP) else 0

  flr <- 1 - precision
  mmr <- fnr

  list(FLR = flr,
       MMR = mmr)

}

#' @title Confusion Matrix for Record Linkage
#'
#' @noRd
get_confusion <- function(TP, FP, FN, TN) {

  cm <- matrix(c(TP, FP, FN, TN), nrow = 2)
  colnames(cm) <- c("Predicted Positive",
                    "Predicted Negative")
  rownames(cm) <- c("Actual Positive",
                    "Actual Negative")
  return(cm)
}
