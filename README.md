
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
# Loading packages
library(fixest)
#> Warning: package 'fixest' was built under R version 4.2.3
library(late.rest)

# Generate a simulated dataset
data <- sim.data()

# Create a predicted compliance score
data2 <- create.score.groups(data, treat = "d", instrument = "z", controls = ~x,
                             n_groups = 10)
#> Warning in doTryCatch(return(expr), name, parentenv, handler): Added noise of
#> the order of 2^-30 to break score ties. Consider using less groups.

# Print summary of predicted compliance scores
summary(data2$pred_p)
#>       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#> -0.0030744  0.0000027  0.0021427  0.2568039  0.7540144  1.0008231

# Run a standard IV
reg1 <- feols(data = data2, y ~ 1 | d ~ z)

# Run basic test-and-select method
reg2 <- run.late.rest(data = data2, yname = "y", treat = "d", instrument = "z", controls = ~g)
#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable
#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

# Run basic test-and-select method
reg3 <- run.late.rest(data = data2, yname = "y", treat = "d", instrument = "z",
                                   controls = ~score_g)
#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

#> Warning in summary.lm(x): essentially perfect fit: summary may be unreliable

# Comparing both regressions
etable(reg1, reg2, reg3)
#>                             reg1             reg2             reg3
#> Dependent Var.:                y                y                y
#>                                                                   
#> Constant        -0.2347 (0.1778) -0.0296 (0.0798) -0.0878 (0.0848)
#> d                0.3507 (0.3369)  0.1318 (0.1400)  0.2090 (0.1327)
#> _______________ ________________ ________________ ________________
#> S.E. type                    IID Heterosked.-rob. Heterosked.-rob.
#> Observations               1,000              300              300
#> R2                      -0.10588          0.00632          0.00399
#> Adj. R2                 -0.10699          0.00299          0.00065
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
