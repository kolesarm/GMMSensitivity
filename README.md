[![Travis build status](https://travis-ci.org/kolesarm/GMMSensitivity.svg?branch=master)](https://travis-ci.org/kolesarm/GMMSensitivity) [![Coverage status](https://codecov.io/gh/kolesarm/GMMSensitivity/branch/master/graph/badge.svg)](https://codecov.io/github/kolesarm/GMMSensitivity?branch=master)

# GMMSensitivity

Estimators and confidence intervals for sensitivity analysis in moment condition
models using procedures from [Armstrong and Kolesár
(2020)](https://arxiv.org/abs/1808.07387).

See the [package vignette](doc/GMMSensitivity.pdf) for a description of the package
(available through `vignette("GMMSensitivity")` once package is installed), and
the package [manual](doc/manual.pdf) for documentation of the package functions.

This software package is based upon work supported by the National Science
Foundation under grant numbers SES-1628939 (Armstrong) and SES-1628878
(Kolesár).

## Installation

You can install the development version from GitHub:
``` r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("kolesarm/GMMSensitivity")
```

- This packages uses the `CVXR` package for testing the package functions (the
  `CVXR` package is not needed to actually use `GMMSensitivity`). If you'd like
  to run these tests locally, you therefore need to install the `CVXR` package.
  In turn, the `CVXR` package package depends on the `Rmpfr` package (it's a
  dependency of the `CVXR` package that this package uses). In order to install
  this package, the system library `libmpfr-dev` needs to be installed. On
  Ubuntu/Debian, the library can be installed by running `sudo apt-get install
  libmpfr-dev`.
