# MRAlasso

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20843253.svg)](https://doi.org/10.5281/zenodo.20843253)

`MRAlasso` is an R package for adaptive penalized Mendelian randomization using summary-level data. The package implements adaptive MR-Lasso (MR-ALasso) and a bootstrap-smoothed MR-ALasso estimator for two-sample Mendelian randomization analyses in the presence of horizontal pleiotropy.


## Installation

You can install the package directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("Qasim-stat/MRAlasso")
```

Then load the package:

```r
library(MRAlasso)
```

## Main function

The main function is

```r
mr_alasso(object)
```

where `object` is an `MRInput` object created using `MendelianRandomization::mr_input()`.

The MR-ALasso estimator can be fitted using

```r
fit = mr_alasso(object)
```

The bootstrap-smoothed version can be fitted using

```r
fit_b = mr_alasso(object, bootstrap = TRUE, B = 200)
```

The archived release is available at

<https://doi.org/10.5281/zenodo.20843253>

