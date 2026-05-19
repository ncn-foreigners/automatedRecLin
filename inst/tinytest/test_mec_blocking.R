library(automatedRecLin)
library(data.table)

options("text2vec.mc.cores" = 1L)

controls_blocking <- list(
  representation = "custom_matrix",
  ann = "kd",
  distance = "euclidean",
  seed = 1
)

expect_blocking_output_contract <- function(result) {
  expect_false(any(c("flr_est", "mmr_est") %in% names(result)))
  expect_false(any(c("training_rule", "training_blocks", "controls_blocking") %in% names(result)))
  expect_equal(names(result$M_est), c("a", "b", "block", "ratio"))
  expect_equal(result$ratio_orientation, "u_over_m")
  expect_equal(result$pooled_model$ratio_orientation, "u_over_m")
  expect_equal(result$n_M_est, NROW(result$M_est))
  expect_equal(result$n_U_est, result$candidate_pair_count - result$n_M_est)
  expect_equal(result$n_U_min, result$candidate_pair_count - result$nu)
  expect_true(result$n_U_est >= result$n_U_min)
  expect_false(any(c("nonmatch_sample_size", "nonmatch_sampling_seed", "prob_ratio") %in% names(result)))
  expect_equal(
    names(result$block_estimates),
    c("block", "n_A", "n_B", "pair_count", "nonmatches_min", "n_M_est", "selected_pairs")
  )
  expect_true(all(result$block_estimates[["n_M_est"]] >= 0))
  expect_true(all(result$block_estimates[["selected_pairs"]] >= 0))
  expect_true(all(result$block_estimates[["n_M_est"]] <=
                    pmin(result$block_estimates[["n_A"]], result$block_estimates[["n_B"]])))
  expect_true(all(result$block_estimates[["selected_pairs"]] <=
                    pmin(result$block_estimates[["n_A"]], result$block_estimates[["n_B"]])))
}

expect_equal(
  automatedRecLin:::select_inverted_mec_indices(
    a = c(1L, 1L, 2L, 2L),
    b = c(1L, 2L, 1L, 2L),
    block = c(1, 1, 1, 1),
    ratio = c(2, 0.1, 0.2, Inf),
    n_M = 2
  ),
  c(2L, 3L)
)
expect_equal(
  automatedRecLin:::select_inverted_mec_indices(
    a = c(1L, 2L),
    b = c(1L, 2L),
    block = c(1, 1),
    ratio = c(NA_real_, 0.5),
    n_M = 1
  ),
  2L
)
expect_equal(
  automatedRecLin:::select_inverted_mec_indices(
    a = c(1L, 2L),
    b = c(1L, 2L),
    block = c(1, 1),
    ratio = c(-1, 0.3),
    n_M = 1
  ),
  2L
)
expect_equal(
  automatedRecLin:::select_inverted_mec_indices(
    a = c(1L, 2L),
    b = c(1L, 2L),
    block = c(1, 1),
    ratio = c(0.1, 0.2),
    n_M = 0
  ),
  integer()
)

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
expect_equal(fit_binary$ratio_orientation, "u_over_m")
expect_equal(fit_binary$n_M_est, 5L)
expect_equal(fit_binary$n_U_est, 0L)
expect_equal(fit_binary$n_U_min, 0L)
expect_equal(fit_binary$nu, 5L)
expect_equal(fit_binary$candidate_pair_count, 5L)
expect_equal(fit_binary$pooled_model$convergence_reason, "structural_no_nonmatch_complement")
expect_equal(
  fit_binary$M_est[, .(a, b, block)],
  data.table(a = 1:5, b = 1:5, block = as.numeric(1:5))
)
expect_true(all(is.finite(fit_binary$M_est[["ratio"]])))
expect_true(all(fit_binary$M_est[["ratio"]] >= 0))
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
expect_blocking_output_contract(fit_binary)
expect_equal(fit_binary$block_estimates[["n_M_est"]], rep(1L, 5L))
expect_equal(fit_binary$block_estimates[["selected_pairs"]], rep(1L, 5L))
expect_null(fit_binary$blocking_result)
expect_null(fit_binary$training_Omega)

fit_binary_print <- paste(capture.output(print(fit_binary)), collapse = "\n")
expect_false(grepl("Estimated false link rate", fit_binary_print))
expect_false(grepl("Estimated missing match rate", fit_binary_print))
expect_true(grepl("Evaluation metrics", fit_binary_print))

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
expect_blocking_output_contract(fit_cpar)

singleton_nonmatch_A <- data.frame(
  name = c("same", "left"),
  surname = c("person", "alpha")
)
singleton_nonmatch_B <- data.frame(
  name = c("same", "right"),
  surname = c("person", "beta")
)
singleton_nonmatch_blocking <- matrix(c(1, 2), ncol = 1)

set.seed(1)
fit_singleton_nonmatch <- mec_blocking(
  A = singleton_nonmatch_A,
  B = singleton_nonmatch_B,
  variables = c("name", "surname"),
  blocking_x = singleton_nonmatch_blocking,
  blocking_y = singleton_nonmatch_blocking,
  controls_blocking = controls_blocking
)

expect_blocking_output_contract(fit_singleton_nonmatch)
expect_equal(fit_singleton_nonmatch$n_M_est, 2L)
expect_equal(fit_singleton_nonmatch$n_U_est, 0L)
expect_equal(fit_singleton_nonmatch$pooled_model$convergence_reason, "structural_no_nonmatch_complement")
expect_equal(fit_singleton_nonmatch$M_est[, .(a, b)], data.table(a = 1:2, b = 1:2))
expect_equal(
  unkey(fit_singleton_nonmatch$block_estimates[, .(block, n_M_est, selected_pairs)]),
  data.table(block = as.numeric(1:2), n_M_est = c(1L, 1L), selected_pairs = c(1L, 1L))
)

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
  true_matches = data.frame(a = 1:6, b = 1:6),
  keep_training_data = TRUE,
  keep_blocking_result = TRUE
)

expect_equal(fit_threshold$ratio_orientation, "u_over_m")
expect_equal(fit_threshold$candidate_pair_count, 6L)
expect_equal(fit_threshold$nu, 2L)
expect_equal(fit_threshold$n_U_min, 4L)
expect_equal(fit_threshold$n_U_est, 4L)
expect_equal(fit_threshold$n_M_est, 2L)
expect_equal(fit_threshold$block_summary[["block"]], as.numeric(1:2))
expect_equal(NROW(fit_threshold$training_Omega), fit_threshold$candidate_pair_count)
expect_true(all(c("gamma_name", "gamma_surname", "init_disagreement", "ratio", "q_est") %in%
                  names(fit_threshold$training_Omega)))
expect_true(all(fit_threshold$training_Omega[["q_est"]] >= 0 &
                  fit_threshold$training_Omega[["q_est"]] <= 1))
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
expect_blocking_output_contract(fit_threshold)

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
