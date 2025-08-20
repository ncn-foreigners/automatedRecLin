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

#' @noRd
get_metrics <- function(TP, FP, FN, TN) {

  recall <- if (TP + FN != 0) TP / (TP + FN) else 0
  precision <- if (TP + FP != 0) TP / (TP + FP) else 0
  fpr <- if (FP + TN != 0) FP / (FP + TN) else 0
  fnr <- if (FN + TP != 0) FN / (FN + TP) else 0
  accuracy <- if (TP + FP + FN + TN != 0) (TP + TN) / (TP + FP + FN + TN) else 0
  specificity <- if (TN + FP != 0) TN / (TN + FP) else 0
  f1_score <- if (precision + recall != 0) 2 * (precision * recall) / (precision + recall) else 0

  list(recall = recall,
       precision = precision,
       fpr = fpr,
       fnr = fnr,
       accuracy = accuracy,
       specificity = specificity,
       f1_score = f1_score)

}

get_confusion <- function(TP, FP, FN, TN) {

  cm <- matrix(c(TP, FP, FN, TN), nrow = 2)
  colnames(cm) <- c("Predicted Positive",
                    "Predicted Negative")
  rownames(cm) <- c("Actual Positive",
                    "Actual Negative")
  return(cm)
}
