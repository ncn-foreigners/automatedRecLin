# Changelog

## Version 1.1.1

- Improved
  [`mec_blocking()`](https://ncn-foreigners.github.io/automatedRecLin/reference/mec_blocking.md)
  by using inverted unsupervised MEC.
- Added `alpha` in
  [`mec_blocking()`](https://ncn-foreigners.github.io/automatedRecLin/reference/mec_blocking.md)
  for controlling the FLR-MMR trade-off.

## Version 1.1.0

CRAN release: 2026-05-08

- Added
  [`mec_blocking()`](https://ncn-foreigners.github.io/automatedRecLin/reference/mec_blocking.md)
  for blocked unsupervised MEC with pooled training and blockwise
  prediction using the `blocking` package.
- Added support for creating comparison vectors on a supplied table of
  record pairs through the `pairs` argument in
  [`comparison_vectors()`](https://ncn-foreigners.github.io/automatedRecLin/reference/comparison_vectors.md).
- Added `census` and `cis` example datasets for larger record linkage
  examples.
- Added a vignette showing MEC with blocking on the `cis` and `census`
  datasets.
- Added optional progress messages via the `verbose` argument in
  [`mec()`](https://ncn-foreigners.github.io/automatedRecLin/reference/mec.md),
  [`train_rec_lin()`](https://ncn-foreigners.github.io/automatedRecLin/reference/train_rec_lin.md),
  [`predict.rec_lin_model()`](https://ncn-foreigners.github.io/automatedRecLin/reference/predict.rec_lin_model.md),
  and
  [`mec_blocking()`](https://ncn-foreigners.github.io/automatedRecLin/reference/mec_blocking.md).
- Improved validation of supplied match and pair tables, including
  clearer checks for row indices, duplicate pairs, missing values, and
  non-finite comparison values.
- Improved print methods for linkage results, including consistent
  percentage formatting for error rates.

## Version 1.0.1

CRAN release: 2025-12-13

- Fixed CRAN errors.

## Version 1.0.0

CRAN release: 2025-11-18

- Implemented comparison functions
  [`abs_distance()`](https://ncn-foreigners.github.io/automatedRecLin/reference/abs_distance.md)
  and
  [`jarowinkler_complement()`](https://ncn-foreigners.github.io/automatedRecLin/reference/jarowinkler_complement.md).
- Added support for comparing two datasets using comparison functions.
- Added support for training a supervised record linkage model using
  probability or density ratio estimation, based on the following
  methods: `"binary"`, `"continuous_parametric"`, and
  `"continuous_nonparametric"`.
- Added support for creating a supervised record linkage model using a
  custom machine learning (ML) classifier.
- Added support for predicting matches based on a record linkage model.
- Added the unsupervised maximum entropy classification (MEC) algorithm
  for record linkage. Supported methods are: `"binary"`,
  `"continuous_parametric"`, `"continuous_nonparametric"`, and
  `"hit_miss"`.
- Added support for creating the predicted set of matches based on: its
  estimated size, a target false link rate (FLR) or a target missing
  match rate (MMR).
- Implemented S3 methods for printing.
- Added support for evaluation when true matches are known.
- Added documentation and examples.
