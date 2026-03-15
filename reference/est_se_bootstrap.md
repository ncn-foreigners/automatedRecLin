# Parametric Bootstrap for Standard Error Estimation in MEC

**WORK IN PROGRESS**

*This function is currently under development.*

## Usage

``` r
est_se_bootstrap(mec_result, B = 100, alpha = 0.05)
```

## Arguments

- mec_result:

  An object of class `mec_rec_lin`.

- B:

  A number of bootstrap iterations.

- alpha:

  A significance level for calculating the confidence interval. Default
  is 0.05 (which yields a 95% confidence interval).

## Author

Adam Struzik
