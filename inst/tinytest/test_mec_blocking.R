library(automatedRecLin)
library(data.table)

options("text2vec.mc.cores" = 1L)

controls_blocking <- list(
  representation = "custom_matrix",
  ann = "kd",
  distance = "euclidean",
  seed = 1
)

expect_rate_bounds <- function(result) {
  expect_true(is.finite(result$flr_est))
  expect_true(result$flr_est >= 0 && result$flr_est <= 1)
  expect_true(is.finite(result$mmr_est))
  expect_true(result$mmr_est >= 0 && result$mmr_est <= 1)
}

unkey <- function(x) {
  data.table::setkey(x, NULL)
  x
}

make_doc_example_data <- function() {
  list(
    A = data.frame(
      name = c("Emma", "Liam", "Olivia", "Noah", "Ava"),
      surname = c("Smith", "Jones", "Brown", "Davis", "Miller"),
      city = c("Boston", "Boston", "Austin", "Austin", "Denver")
    ),
    B = data.frame(
      name = c("Emma", "Liam", "Olivia", "Noah", "Ava"),
      surname = c("Smith", "Jones", "Brown", "Davis", "Miller"),
      city = c("Boston", "Boston", "Austin", "Austin", "Denver")
    ),
    blocking_x = matrix(
      c(1, 0, 0, 1, 1, 1, 2, 0, 0, 2),
      ncol = 2,
      byrow = TRUE
    ),
    true_matches = data.frame(a = 1:5, b = 1:5)
  )
}

confusion_all_links <- matrix(c(5, 0, 0, 20), nrow = 2, ncol = 2)
rownames(confusion_all_links) <- c("Actual Positive", "Actual Negative")
colnames(confusion_all_links) <- c("Predicted Positive", "Predicted Negative")

doc_data <- make_doc_example_data()

set.seed(1)
fit_binary <- mec_blocking(
  A = doc_data$A,
  B = doc_data$B,
  variables = c("name", "surname", "city"),
  blocking_x = doc_data$blocking_x,
  blocking_y = doc_data$blocking_x,
  controls_blocking = controls_blocking,
  true_matches = doc_data$true_matches
)

expect_inherits(fit_binary, "mec_blocking")
expect_equal(fit_binary$training_rule, "all_blocks")
expect_equal(
  fit_binary$M_est[, .(a, b, block)],
  data.table(a = 1:5, b = 1:5, block = as.numeric(1:5))
)
expect_true(all(is.finite(fit_binary$M_est[["ratio"]])))
expect_true(all(fit_binary$M_est[["ratio"]] > 0))
expect_equal(
  unkey(fit_binary$block_summary[, .(block, n_A, n_B, pair_count, nonmatches_min)]),
  data.table(
    block = as.numeric(1:5),
    n_A = rep(1L, 5),
    n_B = rep(1L, 5),
    pair_count = rep(1L, 5),
    nonmatches_min = rep(0L, 5)
  )
)
expect_equal(unlist(fit_binary$block_summary[["A"]]), 1:5)
expect_equal(unlist(fit_binary$block_summary[["B"]]), 1:5)
expect_equal(fit_binary$excluded_records$A, integer())
expect_equal(fit_binary$excluded_records$B, integer())
expect_equal(
  fit_binary$blocking_eval,
  data.table(
    true_matches = 5L,
    preserved_matches = 5L,
    lost_matches = 0L,
    blocking_recall = 1,
    blocking_fnr = 0,
    blocked_pairs = 5L,
    full_pairs = 25L
  )
)
expect_equal(fit_binary$eval_metrics, c(FLR = 0, MMR = 0))
expect_equal(fit_binary$confusion, confusion_all_links)
expect_rate_bounds(fit_binary)
expect_null(fit_binary$blocking_result)
expect_null(fit_binary$training_Omega)

set.seed(1)
fit_cpar <- mec_blocking(
  A = doc_data$A,
  B = doc_data$B,
  variables = c("name", "surname", "city"),
  comparators = list(
    name = jarowinkler_complement(),
    surname = jarowinkler_complement()
  ),
  methods = list(
    name = "continuous_parametric",
    surname = "continuous_parametric"
  ),
  blocking_x = doc_data$blocking_x,
  blocking_y = doc_data$blocking_x,
  controls_blocking = controls_blocking,
  true_matches = doc_data$true_matches
)

expect_equal(fit_cpar$cpar_vars, c("gamma_name", "gamma_surname"))
expect_equal(fit_cpar$b_vars, "gamma_city")
expect_null(fit_cpar$hm_vars)
expect_null(fit_cpar$cnonpar_vars)
expect_null(fit_cpar$ratio_kliep)
expect_null(fit_cpar$ratio_kliep_list)
expect_equal(
  fit_cpar$M_est[, .(a, b, block)],
  data.table(a = 1:5, b = 1:5, block = as.numeric(1:5))
)
expect_equal(fit_cpar$cpar_params[["variable"]], c("gamma_name", "gamma_surname"))
expect_true(all(is.finite(unlist(fit_cpar$cpar_params[, -"variable"]))))
expect_equal(fit_cpar$b_params[["variable"]], "gamma_city")
expect_true(all(is.finite(fit_cpar$b_params[["theta"]])))
expect_true(all(is.finite(fit_cpar$b_params[["eta"]])))
expect_rate_bounds(fit_cpar)

threshold_A <- data.frame(
  name = c("A1", "A2", "A3", "B1", "B2", "B3"),
  surname = c("S1", "S2", "S3", "T1", "T2", "T3")
)
threshold_B <- threshold_A
threshold_blocking <- matrix(c(1, 1, 1, 2, 2, 2), ncol = 1)

confusion_threshold <- matrix(c(2, 0, 4, 30), nrow = 2, ncol = 2)
rownames(confusion_threshold) <- c("Actual Positive", "Actual Negative")
colnames(confusion_threshold) <- c("Predicted Positive", "Predicted Negative")

set.seed(1)
fit_threshold <- mec_blocking(
  A = threshold_A,
  B = threshold_B,
  variables = c("name", "surname"),
  blocking_x = threshold_blocking,
  blocking_y = threshold_blocking,
  controls_blocking = controls_blocking,
  min_training_pairs = 3,
  min_training_nonmatches = 2,
  block_sampling_seed = 1,
  nonmatch_sampling_seed = 1,
  true_matches = data.frame(a = 1:6, b = 1:6),
  keep_training_data = TRUE,
  keep_blocking_result = TRUE
)

expect_equal(fit_threshold$training_rule, "threshold_sampling")
expect_equal(fit_threshold$training_blocks[["block"]], 1)
expect_equal(NROW(fit_threshold$training_Omega), 3L)
expect_true(all(c("gamma_name", "gamma_surname") %in% names(fit_threshold$training_Omega)))
expect_true(!is.null(fit_threshold$blocking_result))
expect_true("result" %in% names(fit_threshold$blocking_result))

included_A <- as.integer(vapply(
  fit_threshold$block_summary[["A"]],
  function(x) x[1L],
  numeric(1)
))

expect_equal(NROW(fit_threshold$block_summary), 2L)
expect_equal(fit_threshold$block_summary[["n_A"]], rep(1L, 2L))
expect_equal(fit_threshold$block_summary[["n_B"]], rep(3L, 2L))
expect_equal(
  fit_threshold$M_est[, .(a, b, block)],
  data.table(
    a = included_A,
    b = included_A,
    block = fit_threshold$block_summary[["block"]]
  )
)
expect_equal(sort(fit_threshold$excluded_records$A), setdiff(seq_len(6L), included_A))
expect_equal(fit_threshold$excluded_records$B, integer())
expect_equal(
  fit_threshold$blocking_eval,
  data.table(
    true_matches = 6L,
    preserved_matches = 2L,
    lost_matches = 4L,
    blocking_recall = 1 / 3,
    blocking_fnr = 2 / 3,
    blocked_pairs = 6L,
    full_pairs = 36L
  )
)
expect_equal(fit_threshold$eval_metrics, c(FLR = 0, MMR = 2 / 3))
expect_equal(fit_threshold$confusion, confusion_threshold)
expect_rate_bounds(fit_threshold)

expect_error(
  mec_blocking(
    A = doc_data$A,
    B = doc_data$B,
    variables = c("name", "surname", "city"),
    methods = list(name = "hit_miss"),
    blocking_x = doc_data$blocking_x,
    blocking_y = doc_data$blocking_x,
    controls_blocking = controls_blocking
  )
)

expect_error(
  mec_blocking(
    A = doc_data$A,
    B = doc_data$B,
    variables = c("name", "surname", "city"),
    blocking_x = doc_data$blocking_x,
    controls_blocking = controls_blocking
  )
)

expect_error(
  mec_blocking(
    A = doc_data$A,
    B = doc_data$B,
    variables = c("name", "surname", "city"),
    blocking_x = doc_data$blocking_x,
    blocking_y = doc_data$blocking_x,
    controls_blocking = c(list(x = doc_data$blocking_x), controls_blocking)
  )
)

expect_error(
  mec_blocking(
    A = doc_data$A,
    B = doc_data$B,
    variables = c("name", "surname", "city"),
    blocking_x = doc_data$blocking_x,
    blocking_y = doc_data$blocking_x,
    controls_blocking = controls_blocking,
    nonmatch_sample_size = 26
  )
)

expect_error(
  mec_blocking(
    A = doc_data$A,
    B = doc_data$B,
    variables = c("name", "surname", "city"),
    blocking_x = doc_data$blocking_x,
    blocking_y = doc_data$blocking_x,
    controls_blocking = controls_blocking,
    min_training_pairs = 1
  )
)
