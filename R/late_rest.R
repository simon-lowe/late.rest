#' Main function
#'
#' @import fixest
#' @import data.table
#' @import lmtest
#' @import sandwich
#' @import readr
#' @importFrom stats lm
#'
#' @param data A data.frame containing the necessary variables to run the model.
#' @param yname Name of the dependent variable
#' @param treat Name of the treatment variable
#' @param instrument Name of the instrument
#' @param controls Controls. Either as a formula, eg ~ x1 + x2, or a vector of strings, eg c("x1", "x2")
#' @param n_folds Number of folds
#' @param p_th Group selection threshold
#' @param joint_est Final estimation method
#' @param csl_est CS method
#'
#' @return fixest object
#' @export
#'
#' @examples
#' # Loading packages
#' library(fixest)
#'
#' # Generate a simulated dataset
#' data <- sim.data()
#'
#' # Run a standard IV
#' reg1 <- feols(data = data, y ~ 1 | d ~ z)
#'
#' # Run the test-and-select method
#' reg2 <- run.late.rest(data = data, yname = "y", treat = "d", instrument = "z", controls = ~g)
#'
#' # Comparing both regressions
#' etable(reg1, reg2)

run.late.rest <- function(data, yname,
                          treat, instrument, controls,
                          n_folds = 2, p_th = 0.05,
                          joint_est = TRUE, csl_est = FALSE){

  # Convert data to data.table
  data.table::setDT(data)

  # Load helper function ----------------------------------------------------
  # source("R/utiliy_functions.R")
  rand_n_list_group <- function(n, n_groups){
    # stopifnot(is.numeric(n) & is.numeric(n_groups), n >= n_groups)

    rem <- n %% n_groups
    mult <- n %/% n_groups
    tmp <- sample(rep(seq(1, n_groups), mult), mult*n_groups, replace = FALSE)
    return(c(tmp, sample(seq(1, n_groups), rem)))
  }

  # Parsing -----------------------------------------------------------------

  if(inherits(controls, "formula")) c_vars <- all.vars(controls[[2]])

  # Checks ------------------------------------------------------------------

  if (!inherits(data, "data.frame")) stop("`late.rest` requires a data.frame like object for analysis.")

  # Check number of observations by group
  tmp <- data[, list(n  = .N), by = mget(c_vars)]
  if(any(tmp$n < 10)) warning("Some groups defined by the covariates have less than 10 observations.")


  # Code --------------------------------------------------------------------
  # Creating variables for the R check
  g = p_val = res = split_x = val = NULL

  # Create the group variable
  data[, g := .GRP, by = mget(c_vars)]

  # Create the split
  data[, split_x := rand_n_list_group(.N, n_folds), by = g]

  # Compute the cross-fitted compliance rates

  f1 <- formula(paste0(treat, " ~ ", instrument))

  if(n_folds > 1){
    tmp <- data.table()
    for(i in 1:n_folds){
      tmp2 <- data[, as.list(lmtest::coeftest(lm(data = .SD[split_x != i], formula = f1), vcov = sandwich::vcovHC, type = "HC3")[2, c(1, 4)]), by = g]
      names(tmp2)[2:3] <- c("pc", "p_val")
      tmp2[, split_x := i]
      tmp <- rbind(tmp, tmp2)
    }

    dat_reg <- data.table::merge.data.table(data, tmp, by = c("g", "split_x"))
  }
  if(n_folds == 1){
    tmp <- data[, as.list(lmtest::coeftest(lm(data = .SD, formula = f1), vcov = sandwich::vcovHC, type = "HC3")[2, c(1, 4)]), by = g]
    names(tmp)[2:3] <- c("pc", "p_val")
    tmp[, split_x := 1]

    dat_reg <- data.table::merge.data.table(data, tmp, by = c("g", "split_x"))
  }

  f2 <- formula(paste0(yname, " ~ 1 | ", treat, " ~ ", instrument))

  if(joint_est == TRUE){
    reg <- fixest::feols(data = dat_reg[p_val <= p_th], f2, vcov = "hetero")
    return(reg)
  }
  if(joint_est == FALSE & csl_est == FALSE){
    coefficient <- NULL
    reg <- fixest::feols(data = dat_reg[p_val <= p_th], f2, split = ~split_x, vcov = "hetero")
    return(setDT(coeftable(reg))[coefficient != "(Intercept)"])
  }
  if(joint_est == FALSE & csl_est == TRUE){
    z_dm <- NULL
    dat_reg[, z_dm := get(instrument) - mean(get(instrument))]
    f3 <- formula(paste0(yname, " ~ pc | ", treat, " ~ z_dm:pc"))
    reg <- fixest::feols(data = dat_reg, f3, vcov = "hetero")
    return(reg)
  }
}

#' Simulate data
#'
#' @import data.table
#' @import stats
#' @importFrom mvtnorm rmvnorm
#' @importFrom gtools quantcut
#'
#' @param nsims Number of simulations
#' @param beta Treatment effect scaling
#' @param alpha Correlation of treatment effect with compliance
#' @param n Sample size
#' @param J Number of groups
#' @param rde Covariance between latency to treat and Y0
#' @param sd Variance of latency to treat
#' @param se Variance of Y0
#' @param s_at Always-taker share
#' @param s_nt Never-taker share
#' @param zeta Heteroskedasticity
#' @param s_eta Quality of prediction of covariates for compliance rate
#'
#' @return Data.frame
#' @export
#'
#' @examples
#' # Create a simulated data.frame
#' data <- sim.data()
sim.data <- function(nsims = 1, beta = 1, alpha = 0.5,
                     n = 1e3, J = 10,
                     rde = 0.3, sd = 1, se = 1,
                     s_at = 0.375, s_nt = 0.375,
                     zeta = 0, s_eta = 0.01){

  # Creating variables for the R check
  d0 = d1 = g = nt = at = pc_tmp = sim_id = x = y = y0 = y1 = z = NULL

  # Draw delta and epsilon
  dist_detau <- mvtnorm::rmvnorm(n*nsims, c(0, 0), matrix(c(sd^2, rde*sd*se, rde*sd*se, se^2), nrow = 2, ncol = 2))

  # Draw z
  dat <- data.table::data.table(z = stats::rbinom(n*nsims, 1, 0.5), sim_id = rep(1:nsims, each = n))

  # Implement DGP
  dat[, d0 := stats::pnorm(dist_detau[ ,1]) < s_at]
  dat[, d1 := stats::pnorm(dist_detau[ ,1]) < 1 - s_nt]
  dat[, c := d1 > d0]
  dat[, nt := d1 == 0]
  dat[, at := d0 == 1]
  dat[, d := d0*(1-z) + d1*z]

  # Create covariates
  dat[, x := dist_detau[ ,1] + stats::rnorm(n*nsims, 0, s_eta)]
  dat[, g := gtools::quantcut(x, q = J, labels = FALSE), by = sim_id]
  # dat[, g := cut(x, breaks = quantile(x, probs = 0:J/J), labels = FALSE), by = sim_id]

  # Create outcomes
  dat[, pc_tmp := mean(c), by = list(sim_id, g)]
  dat[, y0 := (1 + zeta*dist_detau[ ,1])*dist_detau[ ,2]]
  dat[, y1 := beta*(alpha*(pc_tmp - mean(pc_tmp)) + (1-alpha)*stats::rnorm(.N, mean = 0, sd = sd(pc_tmp))) + y0]

  dat[, y := d*y1 + (1-d)*y0]

  return(dat)
}
