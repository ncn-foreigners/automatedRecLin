library(automatedRecLin)
library(data.table)

A_join <- data.table(
  id = c("A1", "A2", "A3"),
  name = c("Ann", "Bob", "Cat"),
  block = c("east", "west", "north")
)
B_join <- data.table(
  id = c("B1", "B2", "B3"),
  name = c("Cat", "Ann", "Dan"),
  block = c("north", "east", "south")
)
links_join <- data.table(
  a = c(1L, 3L),
  b = c(2L, 1L),
  block = c(100L, 200L),
  ratio = c(9.5, 4.25)
)
fit_join <- structure(list(M_est = links_join), class = "mec_blocking")

expected_matches <- data.table(
  a = c(1L, 3L),
  b = c(2L, 1L),
  id.a = c("A1", "A3"),
  name.a = c("Ann", "Cat"),
  block.a = c("east", "north"),
  id.b = c("B2", "B1"),
  name.b = c("Ann", "Cat"),
  block.b = c("east", "north")
)
expected_with_metadata <- data.table(
  a = c(1L, 3L),
  b = c(2L, 1L),
  block = c(100L, 200L),
  ratio = c(9.5, 4.25),
  id.a = c("A1", "A3"),
  name.a = c("Ann", "Cat"),
  block.a = c("east", "north"),
  id.b = c("B2", "B1"),
  name.b = c("Ann", "Cat"),
  block.b = c("east", "north")
)

A_before <- copy(A_join)
B_before <- copy(B_join)
links_before <- copy(links_join)

expect_equal(join_records(fit_join, A_join, B_join), expected_matches)
expect_equal(join_records(as.data.frame(links_join), A_join, B_join), expected_matches)
expect_equal(A_join, A_before)
expect_equal(B_join, B_before)
expect_equal(links_join, links_before)

prediction_links <- links_join[, c("a", "b", "ratio")]
prediction_fit <- structure(list(M_est = prediction_links), class = "rec_lin_predictions")
mec_fit <- structure(list(M_est = prediction_links), class = "mec_rec_lin")
expected_predictions <- expected_matches[, c(
  "a", "b", "id.a", "name.a", "block.a", "id.b", "name.b", "block.b"
), with = FALSE]
expect_equal(join_records(prediction_fit, A_join, B_join), expected_predictions)
expect_equal(join_records(mec_fit, A_join, B_join), expected_predictions)

expect_equal(
  join_records(fit_join, A_join, B_join, keep_from_links = FALSE),
  expected_matches
)
expect_equal(
  join_records(fit_join, A_join, B_join, keep_from_links = TRUE),
  expected_with_metadata
)
expect_equal(
  join_records(as.data.frame(links_join), A_join, B_join, keep_from_links = TRUE),
  expected_with_metadata
)
expected_ratio <- data.table(
  a = c(1L, 3L),
  b = c(2L, 1L),
  ratio = c(9.5, 4.25),
  id.a = c("A1", "A3"),
  name.a = c("Ann", "Cat"),
  block.a = c("east", "north"),
  id.b = c("B2", "B1"),
  name.b = c("Ann", "Cat"),
  block.b = c("east", "north")
)
expect_equal(
  join_records(fit_join, A_join, B_join, keep_from_links = "ratio"),
  expected_ratio
)

expected_custom_suffixes <- copy(expected_matches)
setnames(
  expected_custom_suffixes,
  c("id.a", "name.a", "block.a", "id.b", "name.b", "block.b"),
  c("id_A", "name_A", "block_A", "id_B", "name_B", "block_B")
)
expect_equal(
  join_records(fit_join, A_join, B_join, suffixes = c("_A", "_B")),
  expected_custom_suffixes
)

expected_all_A <- rbind(
  expected_matches,
  data.table(
    a = 2L,
    b = NA_integer_,
    id.a = "A2",
    name.a = "Bob",
    block.a = "west",
    id.b = NA_character_,
    name.b = NA_character_,
    block.b = NA_character_
  )
)
expect_equal(join_records(fit_join, A_join, B_join, all_A = TRUE), expected_all_A)

expected_all_B <- rbind(
  expected_matches,
  data.table(
    a = NA_integer_,
    b = 3L,
    id.a = NA_character_,
    name.a = NA_character_,
    block.a = NA_character_,
    id.b = "B3",
    name.b = "Dan",
    block.b = "south"
  )
)
expect_equal(join_records(fit_join, A_join, B_join, all_B = TRUE), expected_all_B)

expected_all <- rbind(
  expected_all_A,
  data.table(
    a = NA_integer_,
    b = 3L,
    id.a = NA_character_,
    name.a = NA_character_,
    block.a = NA_character_,
    id.b = "B3",
    name.b = "Dan",
    block.b = "south"
  )
)
expect_equal(join_records(fit_join, A_join, B_join, all = TRUE), expected_all)

empty_links <- links_join[0]
expected_empty <- expected_matches[0]
expect_equal(join_records(empty_links, A_join, B_join), expected_empty)

expect_error(join_records(data.table(a = 1L, ratio = 1), A_join, B_join))
expect_error(join_records(data.table(a = c(1L, NA_integer_), b = 1:2), A_join, B_join))
expect_error(join_records(data.table(a = 1.5, b = 1), A_join, B_join))
expect_error(join_records(data.table(a = c(1L, 4L), b = c(1L, 1L)), A_join, B_join))
expect_error(join_records(data.table(a = c(1L, 1L), b = c(1L, 1L)), A_join, B_join))
expect_error(join_records(fit_join, A_join, B_join, suffixes = c(".x", ".x")))
expect_error(join_records(fit_join, A_join, B_join, keep_from_links = NA))
expect_error(join_records(fit_join, A_join, B_join, keep_from_links = "score"))
