# Create Comparison Vectors for Record Linkage

Creates comparison vectors between records in two datasets based on
specified variables and comparison functions.

## Usage

``` r
comparison_vectors(A, B, variables, comparators = NULL, matches = NULL)
```

## Arguments

- A:

  A duplicate-free `data.frame` or `data.table`.

- B:

  A duplicate-free `data.frame` or `data.table`.

- variables:

  A character vector of key variables used to create comparison vectors.

- comparators:

  A named list of functions for comparing pairs of records.

- matches:

  Optional. A `data.frame` or `data.table` indicating known matches.

## Value

Returns a list containing:  

- `Omega` – a `data.table` with comparison vectors between all records
  from both datasets, including optional match information,

- `variables` – a character vector of key variables used for comparison,

- `comparators` – a list of functions used to compare pairs of records,

- `match_prop` – proportion of matches in the smaller dataset.

## Details

Consider two datasets: \\A\\ and \\B\\. For each pair of records \\(a,b)
\in \Omega\\, the function creates a comparison vector
\\\pmb{\gamma}\_{ab} =
(\gamma\_{ab}^1,\gamma\_{ab}^2,\ldots,\gamma\_{ab}^K)'\\ based on
specified \\K\\ variables and comparison functions.

## Note

Each comparison function must return another function, which serves as
the actual comparator.

## Author

Adam Struzik

## Examples

``` r
df_1 <- data.frame(
"name" = c("John", "Emily", "Mark", "Anna", "David"),
"surname" = c("Smith", "Johnson", "Taylor", "Williams", "Brown")
)
df_2 <- data.frame(
  "name" = c("Jon", "Emely", "Marc", "Michael"),
  "surname" = c("Smitth", "Jonson", "Tailor", "Henderson")
)
comparators <- list("name" = jarowinkler_complement(),
                    "surname" = jarowinkler_complement())
matches <- data.frame("a" = 1:3, "b" = 1:3)
result <- comparison_vectors(A = df_1, B = df_2, variables = c("name", "surname"),
                             comparators = comparators, matches = matches)
result
#> Comparison based on the following variables: name, surname.
#> ========================================================
#>        a     b gamma_name gamma_surname match
#>    <int> <int>      <num>         <num> <num>
#> 1:     1     1 0.08333333    0.05555556     1
#> 2:     1     2 1.00000000    1.00000000     0
#> 3:     1     3 1.00000000    0.54444444     0
#> 4:     1     4 0.53571429    1.00000000     0
#> 5:     2     1 1.00000000    1.00000000     0
#> 6:     2     2 0.13333333    0.04761905     1
```
