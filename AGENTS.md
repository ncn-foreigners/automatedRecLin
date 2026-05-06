# AGENTS.md

This file includes guidelines for the Codex app.

## Agent Role

Act as an expert R developer and data scientist
specialized in entity resolution and CRAN package
development.

## Overview of the `automatedRecLin` R package

This is an R package for record linkage based on an
entropy-maximizing classifier. The method is described
in two papers:

- `papers/Lee_et_al_2022.pdf` -- introduces
the method in the binary setting;
- `papers/Struzik_Beresewicz_2026.pdf` -- extends the approach
by including continuous comparisons.

## Code Structure

The core code is included in the `R/` folder, with the
following files:

- `R/bootstrap.R` -- work-in-progress functions for
estimating standard errors of the number of matches estimate;
- `R/comparators.R` -- functions implementing some comparison
methods for records;
- `R/comparison_vectors.R` -- functions for creating comparison
vectors between records in two datasets;
- `R/controls.R` -- controls for the Kullback-Leibler
Importance Procedure used in the package;
- `R/data.R` -- description of datasets included in
the package;
- `R/eval.R` -- functions for evaluating record linkage;
- `R/internals.R` -- minor technical functions;
- `R/methods.R` -- methods for printing objects;
- `R/predict.R` -- `predict()` method for the supervised
version of the algorithm;
- `R/supervised_learning.R` -- functions for training
a supervised MEC record linkage model;
- `R/unsupervised_learning.R` -- the unsupervised MEC
algorithm (the most important feature of the package).

## Code Guidelines

When writing code, stick to the following rules:

- If possible, use the `data.table` R package. However,
take into account that CRAN checks sometimes show problems
with names references, so be careful.
- Don't use the `dplyr` and `tidyr` packages.
- Always ask for adding new dependencies.
- Use variable names that are consistent with the current
names.
- When changing existing code, try to introduce only minimal
and necessary changes.
- Always communicate your changes to the existing code
(via chat).
- Include general comments.
- Implement only methods/approaches that you are directly
asked for or that are present in the papers in the repository.
Don't change the methodology on your own.
- Don't create new directories.
- Don't create a file with the compiled package in the directory.
- If not asked, don't touch the files in `inst/tinytest/`.
- After R CMD check, always remove the check folder and the compiled package
from the directory.

## Documentation Guidelines

When writing documentation, stick to the following rules:

- Use `roxygen2` comments.
- Use American English.
- Use proper technical vocabulary.
- Use proper function/package references according
to CRAN policies.
- To refer to functions from the package, use `[function()]`.
- Refer to functions from other packages as `\link[package]{function()}`. Don't use `::`.
- Use backticks instead of `\code{}`.

## The `blocking` R package

We have also developed the [`blocking`](https://github.com/ncn-foreigners/blocking) R package,
which blocks records for record linkage and deduplication based on approximate nearest neighbor
algorithms (ANNs). Refer to the package when you are asked about blocking in record linkage.
The package is also described in `papers/paper-blocking.pdf`.

