
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Impulse

This package implements the [Chechik &
Koller](https://www.ncbi.nlm.nih.gov/pubmed/19193146) impulse model
using TensorFlow to improve scaleability and allow for the introduction
of priors which improve model interpretability. This model describes
timeseries data using two sigmoidal responses which are sufficient to
capture the dynamics of many biological systems. While this model was
formulated to capture biological dynamics, the model is generally
suitable for any kind of saturation behavior described by half-max
value(s) and assymptote(s).

The core functionality of this package is:

  - simulate timecourse parameters and resulting timecourses
  - fit sigmoid and impulse models to timecourses with or without priors
    on kinetic parameters
  - compare sigmoid and impulse models
  - visualize measurements and parametric fits

## The model

![](man/figures/README-unnamed-chunk-2-1.png)<!-- -->

**sigmoid:** \[
v_{inter}*\frac{1}{1 + \exp(-\beta(t - t_{rise}))}
\]

**implulse:** \[
\frac{1}{1 + \exp(-\beta(t - t_{rise}))} * (v_{final} + (v_{inter} - v_{final}))\frac{1}{1 + \exp(\beta(t - t_{fall}))}
\]

## Impulse functionality

### Fitting Data

The primary functionality in this package is fitting parametric models
to user-supplied timecourses. The vignette: fitting-timecourses.Rmd
simulates time series, fits parameters to the models and carries out
model comparison.

### Formulating priors

The most important contribution of this work is aaplying priors to
impulse models since there are natural constraints on parameter values
which should hold (non-negative rates, non-negative times, rise before
fall). When these constraints are violated, a good fit may occur, but
interpretability of timing and effect sizes will be lost. The vignette:
setting\_priors.Rmd describes how to formulate the priors and can be
used to guide the tuning of parameters for other application.

## Installation

The package is under active development though and the latest set of
features can be obtained by installing from this repository using
`devtools`

``` r
devtools::install_github('calico/impulse')
```
