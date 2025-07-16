#' @title Train a record linkage classifier
#'
#' @export
train_rec_lin <- function(
    A,
    B,
    matches,
    variables,
    comparators = NULL,
    methods = NULL) {

  vectors <- comparison_vectors(A = A,
                                B = B,
                                variables = variables,
                                comparators = comparators,
                                matches = matches)
  gamma <- vectors$gamma
  comparators <- vectors$comparators

  if(!is.null(methods)) {
    stopifnot("`methods` should be a list." =
              is.list(methods))
  }

  missing_variables <- variables[!(variables %in% names(methods))]
  methods[missing_variables] <- "binary"

  if (sum(methods == "binary") > 0) {
    binary_variables <- c(paste0("gamma_", names(which(methods == "binary"))), "match")
    gamma_binary <- gamma[, ..binary_variables]
  }

  if (sum(methods == "continuous_parametric") > 0) {
    continuous_parametric_variables <- c(paste0("gamma_", names(which(methods == "continuous_parametric"))), "match")
    gamma_continuous_parametric <- gamma[, ..continuous_parametric_variables]
  }

  if (sum(methods == "continuous_nonparametric") > 0) {
    continuous_nonparametric_variables <- c(paste0("gamma_", names(which(methods == "continuous_nonparametric"))), "match")
    gamma_continuous_nonparametric <- gamma[, ..continuous_nonparametric_variables]
  }

}
