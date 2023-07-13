
<!-- README.md is generated from README.Rmd. Please edit that file -->

# late.rest

<!-- badges: start -->
<!-- badges: end -->

The goal of late.rest is to implement the approach to estimating LATE
based on a Test-and-Select method as described in [Hazard and Lowe
(2023)](https://www.dropbox.com/s/x7kzrggcy9a614x/LATEPS_Working_Paper.pdf?dl=0)

## Installation

You can install the development version of late.rest from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("simon-lowe/late.rest")
```

## Example

The example below shows basic use with simulated data:

``` r
library(late.rest)
library(fixest)
#> Warning: package 'fixest' was built under R version 4.2.3

# Generate a simulated dataset
data <- sim.data()

# Run a standard IV
reg1 <- feols(data = data, y ~ 1 | d ~ z)

# Run the test-and-select method
reg2 <- run.late.rest(data = data, yname = "y", treat = "d", instrument = "z", controls = ~g)

# Comparing both regressions
etable(reg1, reg2)
#>                              reg1             reg2
#> Dependent Var.:                 y                y
#>                                                   
#> Constant        -0.3063. (0.1706) -0.0512 (0.0916)
#> d                0.6120. (0.3325) 0.3218* (0.1556)
#> _______________ _________________ ________________
#> S.E. type                     IID              IID
#> Observations                1,000              400
#> R2                       -0.22735         -0.01044
#> Adj. R2                  -0.22858         -0.01298
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
