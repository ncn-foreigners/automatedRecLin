library(automatedRecLin)

cmp <- abs_distance()
expect_equal(
  cmp(1, 5),
  4
)

cmp_jw <- jarowinkler_complement()
expect_equal(
  cmp_jw("Smith", "Smitth"),
  0.0555555555555556
)
