# Two-Piece Distributions (R package)

[![R package](https://img.shields.io/badge/language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

`twopiece` is an R package for the **family of two-piece distributions**, a
flexible class of unimodal distributions constructed by joining two rescaled
half-densities of a symmetric baseline distribution at the mode. This simple
construction introduces asymmetry while preserving the tail behaviour and shape
of the chosen baseline, making two-piece distributions a tractable and
interpretable alternative to standard symmetric models.

Given a unimodal symmetric (about 0) probability density function $f$ from the
location-scale family, the two-piece density is:

$$s(x;\, \mu, \sigma_1, \sigma_2, \delta) = \frac{2}{\sigma_1 + \sigma_2}
\begin{cases}
f\left(\dfrac{x - \mu}{\sigma_1};\, \delta\right), & x < \mu, \\
f\left(\dfrac{x - \mu}{\sigma_2};\, \delta\right), & x \geq \mu,
\end{cases}$$

where $\mu$ is the mode, $\sigma_1, \sigma_2 > 0$ are scale parameters
governing the left and right tails respectively, and $\delta$ is an optional
shape parameter of the baseline $f$. Setting $\sigma_1 = \sigma_2$ recovers the
symmetric baseline. Common choices of baseline include the normal, Student-t,
Laplace, and logistic distributions, among others.

## Installation

```r
# install.packages("devtools")
devtools::install_github("FJRubio67/twopiece")
library(twopiece)
```

## Main functions

The package provides density (`d`), distribution (`p`), quantile (`q`), and
random sampling (`r`) functions for two-piece distributions with three or four
parameters:

| Function | Description |
|---|---|
| `dtp3` / `ptp3` / `qtp3` / `rtp3` | Two-piece distributions with 3 parameters ($\mu$, $\sigma_1$, $\sigma_2$) |
| `dtp4` / `ptp4` / `qtp4` / `rtp4` | Two-piece distributions with 4 parameters ($\mu$, $\sigma_1$, $\sigma_2$, $\delta$) |

The baseline distribution is passed as a function argument (e.g. `FUN = dnorm`,
`FUN = dt`), making the package compatible with any symmetric unimodal baseline
available in R.

Three parameterisations of asymmetry are supported via the `param` argument:

| `param` | Parameterisation |
|---|---|
| `"tp"` | Two-piece: separate scale parameters $(\sigma_1, \sigma_2)$ |
| `"eps"` | Epsilon-skew: single scale $\sigma$ and skewness parameter $\gamma \in (-1, 1)$ |
| `"isf"` | Inverse scale factors: single scale $\sigma$ and skewness parameter $\gamma \in (0, \infty)$ |

For full documentation: `?dtp3`, `?dtp4`

## Quick example

```r
library(twopiece)

# Two-piece normal: asymmetric with heavier right tail
x <- seq(-5, 10, length.out = 500)
plot(x, dtp3(x, mu = 0, par1 = 1, par2 = 3, FUN = dnorm, param = "tp"),
     type = "l", ylab = "Density", main = "Two-piece Normal")

# Random sampling
sim <- rtp3(1000, mu = 0, par1 = 1, par2 = 3, FUN = rnorm, param = "tp")
```

## Tutorials

- [The family of two-piece distributions](https://rpubs.com/FJRubio/TPD) — mathematical background and examples
- [twopiece R package](https://rpubs.com/FJRubio/twopiece) — package walkthrough

## Related packages

- [DTP](https://github.com/FJRubio67/DTP) — Double two-piece distributions (R)
- [TPSAS](https://github.com/FJRubio67/TPSAS) — Two-piece sinh-arcsinh distributions (R)
- [FlexBinReg](https://github.com/FJRubio67/FlexBinReg) — Binary regression with two-piece link functions

## License

This package is licensed under the [MIT License](LICENSE).
