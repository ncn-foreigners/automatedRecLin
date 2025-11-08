library(automatedRecLin)
library(data.table)

df_1 <- data.frame(
  "name" = c("James", "Emma", "William", "Olivia", "Thomas",
             "Sophie", "Harry", "Amelia", "George", "Isabella"),
  "surname" = c("Smith", "Johnson", "Brown", "Taylor", "Wilson",
                "Davis", "Clark", "Harris", "Lewis", "Walker")
)
df_2 <- data.frame(
  "name" = c("James", "Ema", "Wimliam", "Olivia", "Charlotte",
             "Henry", "Lucy", "Edward", "Alice", "Jack"),
  "surname" = c("Smith", "Johnson", "Bron", "Tailor", "Moore",
                "Evans", "Hall", "Wright", "Green", "King")
)
comparators <- list("name" = jarowinkler_complement(),
                    "surname" = jarowinkler_complement())
matches <- data.frame("a" = 1:4, "b" = 1:4)

expect_equal(
  head(comparison_vectors(A = df_1, B = df_2, variables = c("name", "surname"), matches = matches)$Omega),
  data.table("a" = rep(1, 6), "b" = 1:6, "gamma_name" = c(1, rep(0, 5)),
             "gamma_surname" = c(1, rep(0, 5)), match = c(1, rep(0, 5)))
)

expect_equal(
  head(comparison_vectors(A = df_1, B = df_2, variables = c("name", "surname"),
                          matches = matches, comparators = comparators)$Omega),
  data.table("a" = rep(1, 6), "b" = 1:6,
             "gamma_name" = c(0, 0.4777777777777778567270, 0.5523809523809524169025,
                              1, 0.5629629629629629983256, 1),
             "gamma_surname" = c(0, 0.5523809523809524169025, 1,
                                 0.5444444444444445085907, 1, 1),
             "match" = c(1, rep(0, 5)))
)
