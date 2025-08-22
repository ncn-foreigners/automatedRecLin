#' @import data.table
#' @import reclin2
#'
#' @title Create Comparison Vectors for Record Linkage
#'
#' @author Adam Struzik
#'
#' @description
#' Creates comparison vectors between records in two datasets based on specified variables
#' and comparison functions.
#'
#' @param A A duplicate-free `data.frame` or `data.table`.
#' @param B A duplicate-free `data.frame` or `data.table`.
#' @param variables A character vector of key variables used to create comparison vectors.
#' @param comparators A named list of functions for comparing pairs of records.
#' @param matches Optional. A `data.frame` or `data.table` indicating known matches.
#'
#' @return
#' Returns a list containing:\cr
#' \itemize{
#' \item{`Omega` -- a `data.table` with comparison vectors between all records from both datasets,
#' including optional match information,}
#' \item{`variables` -- a character vector of key variables used for comparison,}
#' \item{`comparators` -- a list of functions used to compare pairs of records.}
#' }
#'
#' @examples
#' df_1 <- data.frame(
#' "name" = c("John", "Emily", "Mark", "Anna", "David"),
#' "surname" = c("Smith", "Johnson", "Taylor", "Williams", "Brown")
#' )
#' df_2 <- data.frame(
#'   "name" = c("Jon", "Emely", "Marc", "Michael"),
#'   "surname" = c("Smitth", "Jonson", "Tailor", "Henderson")
#' )
#' comparators <- list("name" = reclin2::cmp_jarowinkler(),
#'                     "surname" = reclin2::cmp_jarowinkler())
#' matches <- data.frame("a" = 1:3, "b" = 1:3)
#' result <- comparison_vectors(A = df_1, B = df_2, variables = c("name", "surname"),
#'                              comparators = comparators, matches = matches)
#' result
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
  comparators <- comparators[variables]

  data.table::setDT(A)
  data.table::setDT(B)
  A <- A[, variables, with = FALSE]
  B <- B[, variables, with = FALSE]
  # A[, a := .I]
  # B[, b := .I]
  data.table::set(A, j = "a", value = seq_len(nrow(A)))
  data.table::set(B, j = "b", value = seq_len(nrow(B)))

  # Omega <- data.table::CJ(a = A$a, b = B$b)
  Omega <- data.table::CJ(a = A[["a"]], b = B[["b"]])
  data.table::setkey(Omega, NULL)
  A_values <- A[Omega$a, ]
  B_values <- B[Omega$b, ]
  # A_values[, a := NULL]
  # B_values[, b := NULL]
  data.table::set(A, j = "a", value = NULL)
  data.table::set(B, j = "b", value = NULL)


  gamma_names <- paste0("gamma_", variables)
  gamma_list <- lapply(1:K, function(x) {
    variable <- variables[x]
    return(as.numeric(comparators[[x]](A_values[[variable]], B_values[[variable]])))
  })

  Omega[, (gamma_names) := gamma_list]

  if(!is.null(matches)) {
    data.table::setDT(matches)
    # Omega[, match := as.numeric(paste(a, b) %in% paste(matches$a, matches$b))]
    Omega[, match := as.numeric(paste(.SD[["a"]], .SD[["b"]]) %in% paste(matches[["a"]], matches[["b"]]))]
  }

  structure(
    list(
      Omega = data.table(Omega),
      variables = variables,
      comparators = comparators
    ),
    class = "comparison_vectors"
  )

}
