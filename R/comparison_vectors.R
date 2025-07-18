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
    comparators = NULL,
    matches = NULL) {

  stopifnot("`A` should be a data.frame or a data.table." =
            is.data.frame(A) | is.data.table(A))
  stopifnot("`B` should be a data.frame or a data.table." =
            is.data.frame(B) | is.data.table(B))
  stopifnot("`variables` should be a character vector." =
            is.character(variables))

  stopifnot("Not all variables are present in A." =
            all(variables %in% names(A)))
  stopifnot("Not all variables are present in B." =
            all(variables %in% names(B)))

  if (!is.null(comparators)) {
    stopifnot("`comparators` should be a list." =
              is.list(comparators))
  }

  if (!is.null(matches)) {
    stopifnot("`matches` should be a data.frame or a data.table." =
              is.data.frame(matches) | is.data.table(matches))
    stopifnot("`matches` should consist of two columns: a, b." =
              length(colnames(matches)) == 2,
              all(colnames(matches) == c("a", "b")))
  }

  K <- length(variables)

  missing_variables <- variables[!(variables %in% names(comparators))]
  comparators[missing_variables] <- rep(list(reclin2::cmp_identical()),
                                        length(missing_variables))

  data.table::setDT(A)
  data.table::setDT(B)
  A <- A[, ..variables]
  B <- B[, ..variables]
  A[, a := .I]
  B[, b := .I]

  Omega <- data.table::CJ(a = A$a, b = B$b)
  setkey(Omega, NULL)
  A_values <- A[Omega$a, ]
  B_values <- B[Omega$b, ]
  A_values[, a := NULL]
  B_values[, b := NULL]

  gamma_names <- paste0("gamma_", variables)
  gamma_list <- lapply(1:K, function(x) {
    variable <- variables[x]
    return(as.numeric(comparators[[x]](A_values[[variable]], B_values[[variable]])))
  })

  Omega[, (gamma_names) := gamma_list]

  Omega[, match := as.numeric(paste(a, b) %in% paste(matches$a, matches$b))]

  structure(
    list(
      Omega = data.table(Omega),
      variables = variables,
      comparators = comparators
    ),
    class = "comparison_vectors"
  )

}
