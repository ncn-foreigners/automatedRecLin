# Join Records Using Linkage Results

Joins two datasets using row-index pairs returned by record linkage.

## Usage

``` r
join_records(
  links,
  A,
  B,
  all = FALSE,
  all_A = all,
  all_B = all,
  suffixes = c(".a", ".b"),
  keep_from_links = FALSE
)
```

## Arguments

- links:

  A linkage result from
  [`mec()`](https://ncn-foreigners.github.io/automatedRecLin/reference/mec.md),
  [`predict.rec_lin_model()`](https://ncn-foreigners.github.io/automatedRecLin/reference/predict.rec_lin_model.md),
  [`mec_blocking()`](https://ncn-foreigners.github.io/automatedRecLin/reference/mec_blocking.md),
  or a `data.frame`/`data.table` with columns `a` and `b`.

- A:

  A `data.frame` or `data.table`.

- B:

  A `data.frame` or `data.table`.

- all:

  Logical indicating whether to include unmatched records from both
  datasets.

- all_A:

  Logical indicating whether to include unmatched records from `A`.

- all_B:

  Logical indicating whether to include unmatched records from `B`.

- suffixes:

  A character vector of length two used to distinguish columns from `A`
  and `B` when their names conflict with each other or with linkage
  columns.

- keep_from_links:

  Logical or character vector. If `FALSE`, only columns `a` and `b` are
  kept from the linkage table. If `TRUE`, all non-index linkage columns
  are kept. If a character vector, only the selected non-index linkage
  columns are kept.

## Value

Returns a `data.table` containing:  

- `a` – row indices of records from `A`,

- `b` – row indices of records from `B`,

- columns selected from `links` – linkage metadata kept according to
  `keep_from_links`, if requested,

- columns from `A` – values of records from `A`, with `suffixes[1]`
  added when needed,

- columns from `B` – values of records from `B`, with `suffixes[2]`
  added when needed.

## Note

The function follows the general design of
[link()](https://rdrr.io/pkg/reclin2/man/link.html), adjusted to linkage
results used in `automatedRecLin`.

## Author

Adam Struzik

## Examples

``` r
A <- data.frame(name = c("James", "Emma"), age = c(30, 28))
B <- data.frame(name = c("James", "Emily"), city = c("Boston", "Denver"))
links <- data.frame(a = 1, b = 1, ratio = 10)
join_records(links, A, B)
#>        a     b name.a   age name.b   city
#>    <int> <int> <char> <num> <char> <char>
#> 1:     1     1  James    30  James Boston
```
