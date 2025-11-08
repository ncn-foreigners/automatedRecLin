# Changelog

## Version 1.0.0

- Implemented comparison functions `abs_distance` and
  `jarowinkler_complement`.
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
  `"hit-miss"`.
- Added support for creating the predicted set of matches based on: its
  estimated size, a target False Link Rate (FLR) or a target Missing
  Match Rate (MMR).
- Implemented S3 methods for printing.
- Added support for evaluation when true matches are known.
- Added documentation and examples.
