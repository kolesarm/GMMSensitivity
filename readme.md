[![Travis build status](https://travis-ci.org/kolesarm/GMMSensitivity.svg?branch=master)](https://travis-ci.org/kolesarm/GMMSensitivity) [![Coverage status](https://codecov.io/gh/kolesarm/GMMSensitivity/branch/master/graph/badge.svg)](https://codecov.io/github/kolesarm/GMMSensitivity?branch=master)

# GMMSensitivity

Estimators and confidence intervals for sensitivity analysis in moment condition
models using procedures from [Armstrong and Kolesár
(2018)](https://arxiv.org/abs/1808.07387).

See vignette `GMMSensitivityExample.pdf` (stored under `doc/`) for description of
the package.

This software package is based upon work supported by the National Science
Foundation under grant numbers SES-1628939 (Armstrong) and SES-1628878
(Kolesár).

## Installation

You can install the package manually by downloading the source code here, or
using the function `install_github()` from the `devtools` package:

```
install.packages("devtools") ## if devtools package not installed
devtools::install_github("kolesarm/GMMSensitivity")
```

Note: to install the `Rmpfr` package (a dependency of `CVXR` package that this
  package uses), one needs the `libmpfr-dev` library. On Ubuntu/Debian, the
  library can be installed by running `sudo apt-get install libmpfr-dev`
