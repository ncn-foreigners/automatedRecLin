#' @import data.table
#' @importFrom reclin2 cmp_identical
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
#' @details
#' Consider two datasets: \eqn{A} and \eqn{B}.
#' For each pair of records \eqn{(a,b) \in \Omega},
#' the function creates a comparison vector
#' \eqn{\pmb{\gamma}_{ab} = (\gamma_{ab}^1,\gamma_{ab}^2,\ldots,\gamma_{ab}^K)'}
#' based on specified \eqn{K} variables and comparison functions.
#'
#' @note
#' Each comparison function must return another function,
#' which serves as the actual comparator.
#'
#' @return
#' Returns a list containing:\cr
#' \itemize{
#' \item{`Omega` -- a `data.table` with comparison vectors between all records from both datasets,
#' including optional match information,}
#' \item{`variables` -- a character vector of key variables used for comparison,}
#' \item{`comparators` -- a list of functions used to compare pairs of records,}
#' \item{`match_prop` -- proportion of matches in the smaller dataset.}
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
#' comparators <- list("name" = jarowinkler_complement(),
#'                     "surname" = jarowinkler_complement())
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
  stopifnot("`variables` should contain at least one variable." =
              length(variables) > 0L)
  stopifnot("Not all variables are present in A." =
              all(variables %in% names(A)))
  stopifnot("Not all variables are present in B." =
              all(variables %in% names(B)))

  if (!is.null(comparators)) {
    stopifnot("`comparators` should be a list." =
                is.list(comparators))
  }

  K <- length(variables)

  missing_variables <- variables[!(variables %in% names(comparators))]
  comparators[missing_variables] <- rep(list(reclin2::cmp_identical()),
                                        length(missing_variables))
  comparators <- comparators[variables]

  data.table::setDT(A)
  data.table::setDT(B)
  n_A <- nrow(A)
  n_B <- nrow(B)

  if (!is.null(matches)) {
    validate_match_pairs(matches, n_A, n_B)
  }

  Omega <- data.table::CJ(a = seq_len(n_A), b = seq_len(n_B))
  data.table::setkey(Omega, NULL)

  gamma_names <- paste0("gamma_", variables)
  omega_a <- Omega[["a"]]
  omega_b <- Omega[["b"]]

  # Compute one comparison column at a time to avoid copying the full key-variable blocks.
  gamma_list <- lapply(1:K, function(x) {
    variable <- variables[x]
    as.numeric(comparators[[x]](A[[variable]][omega_a], B[[variable]][omega_b]))
  })

  invalid_counts <- vapply(gamma_list, function(gamma) {
    sum(!is.finite(gamma))
  }, numeric(1))

  if (any(invalid_counts > 0)) {
    invalid_vars <- variables[invalid_counts > 0]
    invalid_details <- sprintf(
      "%s (%d invalid value%s)",
      invalid_vars,
      invalid_counts[invalid_counts > 0],
      ifelse(invalid_counts[invalid_counts > 0] == 1, "", "s")
    )
    stop(
      sprintf(
        "Comparison variables produced missing or non-finite values: %s. Please handle missing key values or adjust comparators before running record linkage.",
        paste(invalid_details, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  Omega[, (gamma_names) := gamma_list]

  if (!is.null(matches)) {
    data.table::setDT(matches)
    match_idx <- data.table(
      omega_idx = seq_len(NROW(Omega)),
      a = Omega[["a"]],
      b = Omega[["b"]]
    )[matches[, c("a", "b"), with = FALSE], on = c("a", "b"), nomatch = 0L][["omega_idx"]]
    Omega[, match := 0]
    Omega[match_idx, match := 1]
  }

  structure(
    list(
      Omega = data.table(Omega),
      variables = variables,
      comparators = comparators,
      match_prop = if (is.null(matches)) NULL else NROW(matches) / NROW(Omega) * max(n_A, n_B)
    ),
    class = "comparison_vectors"
  )

}
