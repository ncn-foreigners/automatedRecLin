#' @import data.table
#' @importFrom reclin2 cmp_identical
#'
#' @title Create comparison vectors
#'
#' @export
comparison_vectors <- function(
    A,
    B,
    variables,
    comparators = NULL) {

  stopifnot("Not all variables are present in A." = all(variables %in% names(A)))
  stopifnot("Not all variables are present in B." = all(variables %in% names(B)))

  K <- length(variables)

  missing_variables <- variables[!(variables %in% names(comparators))]
  comparators[missing_variables] <- rep(list(reclin2::cmp_identical()), length(missing_variables))

  A <- A[, ..variables]
  B <- B[, ..variables]
  A[, a := .I]
  B[, b := .I]

  gamma <- data.table::CJ(a = A$a, b = B$b)
  setkey(gamma, NULL)
  A_values <- A[gamma$a, ]
  B_values <- B[gamma$b, ]
  A_values[, a := NULL]
  B_values[, b := NULL]

  gamma_names <- paste0("gamma_", 1:K)
  gamma_list <- lapply(1:K, function(x) {
    variable <- names(comparators)[x]
    return(as.numeric(comparators[[x]](A_values[[variable]], B_values[[variable]])))
  })

  gamma[, (gamma_names) := gamma_list]

  data.table(gamma)

}
