library(data.table)
library(automatedRecLin)
library(blocking)

options(text2vec.mc.cores = 1L)

set.seed(123)
df <- fread("internal/data-new.csv.gz", na.strings = c("<NA>"))
df1 <- df[!duplicated(true_id)]
df2 <- df[duplicated(true_id)]
df1_un <- df1[!(true_id %in% df2[["true_id"]])]
df1 <- df1[true_id %in% df2[["true_id"]]]
to_move <- sample(1:NROW(df1_un), size = 1000)
df2 <- rbind(df2, df1_un[to_move])
df1 <- rbind(df1, df1_un[!to_move])

df1[is.na(df1)] <- ""
df2[is.na(df2)] <- ""

true_matches <- merge(
  x = df1[, .(a = .I, true_id)],
  y = df2[, .(b = .I, true_id)],
  by = "true_id"
)[, .(a, b)]

variables <- c(
  "fname",
  "surname",
  "dob_year",
  "dob_mon",
  "dob_day"#,
  # "country",
  # "date"
)

comparators <- list(
  "fname" = jarowinkler_complement(),
  "surname" = jarowinkler_complement()
)

methods <- list(
  "fname" = "continuous_parametric",
  "surname" = "continuous_parametric"
)

df1[, txt := paste(fname, surname, gsub("/", " ", date), country)]
df2[, txt := paste(fname, surname, gsub("/", " ", date), country)]
df1[, dob_year := substr(date, 1, 4)]
df2[, dob_year := substr(date, 1, 4)]
df1[, dob_mon := substr(date, 6, 7)]
df2[, dob_mon := substr(date, 6, 7)]
df1[, dob_day := substr(date, 9, 10)]
df2[, dob_day := substr(date, 9, 10)]

ann_control_pars <- controls_ann()
ann_control_pars$nnd$epsilon <- 0.5
set.seed(1)
result <- mec_blocking(
  A = df1,
  B = df2,
  variables = variables,
  comparators = comparators,
  methods = methods,
  blocking_x = df1[["txt"]],
  blocking_y = df2[["txt"]],
  true_matches = true_matches,
  verbose = TRUE,
  rho = 1,
  keep_blocking_result = TRUE,
  controls_blocking = list(
    control_ann = ann_control_pars,
    control_txt = controls_txt(n_shingles = 3L)
  )
)

decorate_pair_keys <- function(pairs, A, B, key_vars) {
  pairs <- copy(as.data.table(pairs))

  A_keys <- copy(A[, key_vars, with = FALSE])
  A_keys[, a := .I]
  setcolorder(A_keys, "a")
  setnames(A_keys, key_vars, paste0(key_vars, "_A"))

  B_keys <- copy(B[, key_vars, with = FALSE])
  B_keys[, b := .I]
  setcolorder(B_keys, "b")
  setnames(B_keys, key_vars, paste0(key_vars, "_B"))

  pairs <- A_keys[pairs, on = "a"]
  pairs <- B_keys[pairs, on = "b"]

  key_cols <- as.vector(rbind(paste0(key_vars, "_A"), paste0(key_vars, "_B")))
  leading_cols <- intersect(
    c(
      "a", "b", "block", "preserved_by_blocking",
      "block_n_A", "block_n_B", "block_pair_count",
      "density_ratio", "match_prob_est"
    ),
    names(pairs)
  )
  setcolorder(pairs, c(leading_cols, key_cols, setdiff(names(pairs), c(leading_cols, key_cols))))
  pairs[]
}

make_record_block_maps <- function(block_summary) {
  block_summary <- as.data.table(block_summary)

  if (NROW(block_summary) == 0L) {
    return(list(
      A = data.table(a = integer(), block_A = numeric()),
      B = data.table(b = integer(), block_B = numeric())
    ))
  }

  list(
    A = rbindlist(lapply(seq_len(NROW(block_summary)), function(i) {
      data.table(a = block_summary[["A"]][[i]], block_A = block_summary[["block"]][i])
    })),
    B = rbindlist(lapply(seq_len(NROW(block_summary)), function(i) {
      data.table(b = block_summary[["B"]][[i]], block_B = block_summary[["block"]][i])
    }))
  )
}

add_match_probability <- function(false_links, prob_ratio, pooled_model) {
  false_links[, match_prob_est := NA_real_]

  if (NROW(false_links) == 0L) {
    return(false_links[])
  }

  if (identical(prob_ratio, "1") && !is.null(pooled_model[["prob_est"]])) {
    false_links[, match_prob_est := pmin(pooled_model[["prob_est"]] * density_ratio, 1)]
  } else if (all(c("block_n_M_est", "block_pair_count") %in% names(false_links))) {
    false_links[, match_prob_est := pmin(
      block_n_M_est * density_ratio /
        (block_n_M_est * (density_ratio - 1) + block_pair_count),
      1
    )]
  }

  false_links[!is.finite(match_prob_est), match_prob_est := NA_real_]
  false_links[]
}

pred_matches <- copy(as.data.table(result$M_est))
true_matches_dt <- copy(as.data.table(true_matches))

block_diagnostics <- result$block_estimates[
  ,
  .(
    block,
    block_n_A = n_A,
    block_n_B = n_B,
    block_pair_count = pair_count,
    block_n_M_est = n_M_est
  )
]

false_links <- pred_matches[!true_matches_dt, on = .(a, b)]
setnames(false_links, "ratio", "density_ratio")
false_links <- block_diagnostics[false_links, on = "block"]
false_links <- add_match_probability(false_links, result$prob_ratio, result$pooled_model)
false_links[, block_n_M_est := NULL]
false_links_table <- decorate_pair_keys(false_links, df1, df2, variables)

missed_matches <- true_matches_dt[!pred_matches, on = .(a, b)]
record_block_maps <- make_record_block_maps(result$block_summary)
missed_matches <- record_block_maps$A[missed_matches, on = "a"]
missed_matches <- record_block_maps$B[missed_matches, on = "b"]
missed_matches[, preserved_by_blocking := !is.na(block_A) & !is.na(block_B) & block_A == block_B]
missed_matches[, block := block_A]
missed_matches[preserved_by_blocking == FALSE, block := NA]
missed_matches[, c("block_A", "block_B") := NULL]
missed_matches_table <- decorate_pair_keys(missed_matches, df1, df2, variables)

cat(sprintf("False links: %s\n", format(NROW(false_links_table), big.mark = ",")))
cat(sprintf("Missed matches: %s\n", format(NROW(missed_matches_table), big.mark = ",")))

if (!is.null(result$confusion)) {
  stopifnot(NROW(false_links_table) == result$confusion["Actual Negative", "Predicted Positive"])
  stopifnot(NROW(missed_matches_table) == result$confusion["Actual Positive", "Predicted Negative"])
}

stopifnot(NROW(false_links_table[true_matches_dt, on = .(a, b), nomatch = 0L]) == 0L)
stopifnot(NROW(missed_matches_table[!true_matches_dt, on = .(a, b)]) == 0L)

if (NROW(false_links_table) > 0L) {
  stopifnot(!anyNA(false_links_table[, .(block_n_A, block_n_B, block_pair_count, density_ratio)]))
  stopifnot(all(
    is.na(false_links_table[["match_prob_est"]]) |
      (false_links_table[["match_prob_est"]] >= 0 & false_links_table[["match_prob_est"]] <= 1)
  ))
}
