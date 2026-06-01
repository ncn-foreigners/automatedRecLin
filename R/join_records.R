#' @import data.table
#'
#' @title Join Records Using Linkage Results
#'
#' @author Adam Struzik
#'
#' @description
#' Joins two datasets using row-index pairs returned by record linkage.
#'
#' @note
#' The function follows the general design of \link[reclin2:link]{link()}, adjusted
#' to linkage results used in `automatedRecLin`.
#'
#' @param links A linkage result from [mec()], [predict.rec_lin_model()],
#' [mec_blocking()], or a `data.frame`/`data.table` with columns `a` and `b`.
#' @param A A `data.frame` or `data.table`.
#' @param B A `data.frame` or `data.table`.
#' @param all Logical indicating whether to include unmatched records from both
#' datasets.
#' @param all_A Logical indicating whether to include unmatched records from `A`.
#' @param all_B Logical indicating whether to include unmatched records from `B`.
#' @param suffixes A character vector of length two used to distinguish columns
#' from `A` and `B` when their names conflict with each other or with linkage
#' columns.
#' @param keep_from_links Logical or character vector. If `FALSE`, only columns
#' `a` and `b` are kept from the linkage table. If `TRUE`, all non-index linkage
#' columns are kept. If a character vector, only the selected non-index linkage
#' columns are kept.
#'
#' @return
#' Returns a `data.table` containing:\cr
#' \itemize{
#' \item `a` -- row indices of records from `A`,
#' \item `b` -- row indices of records from `B`,
#' \item columns selected from `links` -- linkage metadata kept according to
#' `keep_from_links`, if requested,
#' \item columns from `A` -- values of records from `A`, with `suffixes[1]`
#' added when needed,
#' \item columns from `B` -- values of records from `B`, with `suffixes[2]`
#' added when needed.
#' }
#'
#' @examples
#' A <- data.frame(name = c("James", "Emma"), age = c(30, 28))
#' B <- data.frame(name = c("James", "Emily"), city = c("Boston", "Denver"))
#' links <- data.frame(a = 1, b = 1, ratio = 10)
#' join_records(links, A, B)
#' @export
join_records <- function(links,
                         A,
                         B,
                         all = FALSE,
                         all_A = all,
                         all_B = all,
                         suffixes = c(".a", ".b"),
                         keep_from_links = FALSE) {

  if (!(is.data.frame(A) || is.data.table(A))) {
    stop("`A` should be a data.frame or a data.table.")
  }
  if (!(is.data.frame(B) || is.data.table(B))) {
    stop("`B` should be a data.frame or a data.table.")
  }
  if (anyDuplicated(names(A)) > 0L) {
    stop("`A` should have unique column names.")
  }
  if (anyDuplicated(names(B)) > 0L) {
    stop("`B` should have unique column names.")
  }

  validate_logical_scalar <- function(x, arg_name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(sprintf("`%s` should be TRUE or FALSE.", arg_name))
    }
  }
  validate_logical_scalar(all, "all")
  validate_logical_scalar(all_A, "all_A")
  validate_logical_scalar(all_B, "all_B")

  if (!is.character(suffixes) || length(suffixes) != 2L ||
      anyNA(suffixes) || any(!nzchar(suffixes)) ||
      suffixes[1] == suffixes[2]) {
    stop("`suffixes` should contain two distinct non-empty character values.")
  }

  if (inherits(links, c("mec_rec_lin", "rec_lin_predictions", "mec_blocking"))) {
    if (is.null(links$M_est)) {
      stop("`links` should contain an `M_est` element.")
    }
    link_table <- links$M_est
  } else if (is.data.frame(links) || is.data.table(links)) {
    link_table <- links
  } else {
    stop("`links` should be a record linkage result or a data.frame/data.table.")
  }

  A <- data.table::copy(data.table::as.data.table(A))
  B <- data.table::copy(data.table::as.data.table(B))
  link_table <- data.table::copy(data.table::as.data.table(link_table))

  if (anyDuplicated(names(link_table)) > 0L) {
    stop("`links` should have unique column names.")
  }
  if (!all(c("a", "b") %in% names(link_table))) {
    stop("`links` should contain columns: a, b.")
  }
  if (anyNA(link_table[["a"]]) || anyNA(link_table[["b"]])) {
    stop("`links` cannot contain missing values in columns a and b.")
  }
  if (!is.numeric(link_table[["a"]]) || !is.numeric(link_table[["b"]])) {
    stop("`links` should contain numeric row indices in columns a and b.")
  }
  if (any(link_table[["a"]] != as.integer(link_table[["a"]])) ||
      any(link_table[["b"]] != as.integer(link_table[["b"]]))) {
    stop("`links` should contain integer row indices in columns a and b.")
  }
  if (any(link_table[["a"]] < 1) || any(link_table[["b"]] < 1)) {
    stop("`links` should contain positive row indices in columns a and b.")
  }
  if (any(link_table[["a"]] > nrow(A)) || any(link_table[["b"]] > nrow(B))) {
    stop("`links` contains row indices outside the input datasets.")
  }
  if (anyDuplicated(link_table[, c("a", "b"), with = FALSE]) > 0L) {
    stop("`links` should not contain duplicate record pairs.")
  }

  data.table::set(link_table, j = "a", value = as.integer(link_table[["a"]]))
  data.table::set(link_table, j = "b", value = as.integer(link_table[["b"]]))

  link_columns <- setdiff(names(link_table), c("a", "b"))
  if (is.logical(keep_from_links)) {
    if (length(keep_from_links) != 1L || is.na(keep_from_links)) {
      stop("`keep_from_links` should be TRUE, FALSE, or a character vector.")
    }
    keep_from_links <- if (keep_from_links) link_columns else character()
  } else if (!is.character(keep_from_links) || anyNA(keep_from_links)) {
    stop("`keep_from_links` should be TRUE, FALSE, or a character vector.")
  } else if (anyDuplicated(keep_from_links) > 0L) {
    stop("`keep_from_links` should not contain duplicate column names.")
  } else if (!all(keep_from_links %in% link_columns)) {
    stop("`keep_from_links` contains columns not present in `links`.")
  }

  link_out <- link_table[, c("a", "b", keep_from_links), with = FALSE]

  A_names <- names(A)
  B_names <- names(B)
  reserved_names <- names(link_out)
  shared_names <- intersect(A_names, B_names)
  A_out_names <- A_names
  B_out_names <- B_names
  A_out_names[A_names %in% c(shared_names, reserved_names)] <-
    paste0(A_out_names[A_names %in% c(shared_names, reserved_names)], suffixes[1])
  B_out_names[B_names %in% c(shared_names, reserved_names)] <-
    paste0(B_out_names[B_names %in% c(shared_names, reserved_names)], suffixes[2])

  output_names <- c(reserved_names, A_out_names, B_out_names)
  if (anyDuplicated(output_names) > 0L) {
    stop("Output column names are not unique. Use different `suffixes` or rename input columns.")
  }

  make_missing_rows <- function(x, n) {
    if (n == 0L) {
      return(x[integer()])
    }
    data.table::setkey(x[rep(NA_integer_, n)], NULL)
  }

  A_linked <- A[link_table[["a"]]]
  B_linked <- B[link_table[["b"]]]
  data.table::setnames(A_linked, A_out_names)
  data.table::setnames(B_linked, B_out_names)

  result <- cbind(link_out, A_linked, B_linked)

  extra_rows <- list()
  if (isTRUE(all_A)) {
    unmatched_A <- setdiff(seq_len(nrow(A)), link_table[["a"]])
    link_A <- make_missing_rows(link_out, length(unmatched_A))
    data.table::set(link_A, j = "a", value = unmatched_A)
    A_extra <- A[unmatched_A]
    B_extra <- make_missing_rows(B, length(unmatched_A))
    data.table::setnames(A_extra, A_out_names)
    data.table::setnames(B_extra, B_out_names)
    extra_rows[[length(extra_rows) + 1L]] <- cbind(link_A, A_extra, B_extra)
  }
  if (isTRUE(all_B)) {
    unmatched_B <- setdiff(seq_len(nrow(B)), link_table[["b"]])
    link_B <- make_missing_rows(link_out, length(unmatched_B))
    data.table::set(link_B, j = "b", value = unmatched_B)
    A_extra <- make_missing_rows(A, length(unmatched_B))
    B_extra <- B[unmatched_B]
    data.table::setnames(A_extra, A_out_names)
    data.table::setnames(B_extra, B_out_names)
    extra_rows[[length(extra_rows) + 1L]] <- cbind(link_B, A_extra, B_extra)
  }

  if (length(extra_rows) > 0L) {
    result <- data.table::rbindlist(c(list(result), extra_rows), use.names = TRUE)
  }

  result[]
}
