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
#' @param one_sided_test Use of a one-sided t-test for selection.
#' @param p_th Group selection threshold
#' @param joint_est Final estimation method
#' @param csl_est CS method
#' @param inter_drop "Dropping" through interactions
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
#' reg2 <- run.late.rest(dat = data, yname = "y", treat = "d", instrument = "z", controls = ~g)
#'
#' # Comparing both regressions
#' etable(reg1, reg2)

run.late.rest <- function(data, yname,
                          treat, instrument, controls,
                          n_folds = 2, one_sided_test = TRUE, p_th = 0.05,
                          joint_est = TRUE, csl_est = FALSE,
                          inter_drop = FALSE){

  # Convert data to data.table
  # data.table::setDT(data)
  dat <- as.data.table(data)

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

  if (!inherits(dat, "data.frame")) stop("`late.rest` requires a data.frame like object for analysis.")

  # Check number of observations by group
  tmp <- dat[, list(n  = .N), by = mget(c_vars)]
  if(any(tmp$n < 10)) warning("Some groups defined by the covariates have less than 10 observations.")


  # Code --------------------------------------------------------------------
  # Creating variables for the R check
  g = p_val = t_val = res = split_x = val = NULL
  keep_g = NULL

  # Create the group variable
  dat[, g := .GRP, by = mget(c_vars)]

  # Create the split
  repeat_split <- TRUE
  while(repeat_split == TRUE){
    dat[, split_x := rand_n_list_group(.N, n_folds), by = g]

    z_by_g <- NULL
    test_split <- dat[, list(z_by_g = mean(get(instrument))), by = c("g", "split_x")]
    if(any(test_split$z_by_g %in% c(0,1))){
      warning("Had to repeat split because of perfect imbalance.")
    } else {
      repeat_split <- FALSE
    }
  }


  # Compute the cross-fitted compliance rates

  f1 <- formula(paste0(treat, " ~ ", instrument))

  if(n_folds > 1){
    tmp <- data.table()
    for(i in 1:n_folds){
      tmp2 <- dat[, as.list(lmtest::coeftest(lm(data = .SD[split_x != i], formula = f1), vcov = sandwich::vcovHC, type = "HC3")[2, c(1, 3)]), by = g]
      names(tmp2)[2:3] <- c("pc", "t_val")
      tmp2[, split_x := i]
      tmp <- rbind(tmp, tmp2)
      # tmp[is.nan(t_val), t_val := 0]
    }

    dat_reg <- data.table::merge.data.table(dat, tmp, by = c("g", "split_x"))
  }
  if(n_folds == 1){
    tmp <- dat[, as.list(lmtest::coeftest(lm(data = .SD, formula = f1), vcov = sandwich::vcovHC, type = "HC3")[2, c(1, 3)]), by = g]
    names(tmp)[2:3] <- c("pc", "t_val")
    tmp[, split_x := 1]

    dat_reg <- data.table::merge.data.table(dat, tmp, by = c("g", "split_x"))
  }

  f2 <- formula(paste0(yname, " ~ 1 | ", treat, " ~ ", instrument))

  if(joint_est == TRUE){
    if(one_sided_test == TRUE){
      reg <- fixest::feols(data = dat_reg[t_val >= qnorm(1-p_th)], f2, vcov = "hetero")
    }
    if(one_sided_test == FALSE){
      reg <- fixest::feols(data = dat_reg[abs(t_val) >= qnorm(1-p_th/2)], f2, vcov = "hetero")
    }
    return(reg)
  }
  # if(joint_est == FALSE & csl_est == FALSE){
  if(joint_est == FALSE & csl_est == FALSE & inter_drop == FALSE){
    coefficient <- NULL
    # reg <- fixest::feols(data = dat_reg[p_val <= p_th], f2, split = ~split_x, vcov = "hetero")
    if(one_sided_test == TRUE){
      reg <- fixest::feols(data = dat_reg[t_val >= qnorm(1-p_th)], f2, split = ~split_x, vcov = "hetero")
    }
    if(one_sided_test == FALSE){
      reg <- fixest::feols(data = dat_reg[abs(t_val) >= qnorm(1-p_th/2)], f2, split = ~split_x, vcov = "hetero")
    }
    return(setDT(coeftable(reg))[coefficient != "(Intercept)"])
  }
  if(joint_est == FALSE & csl_est == TRUE){
    z_dm <- NULL
    pc <- NULL
    int_z_dm_pc <- NULL
    dat_reg[, z_dm := get(instrument) - mean(get(instrument))]
    dat_reg[, int_z_dm_pc := z_dm * pc]
    # f3 <- formula(paste0(yname, " ~ pc | ", treat, " ~ z_dm:pc"))
    f3 <- formula(paste0(yname, " ~ pc | ", treat, " ~ int_z_dm_pc"))
    reg <- fixest::feols(data = dat_reg, f3, vcov = "hetero")
    return(reg)
  }
  if(joint_est == FALSE & inter_drop == TRUE){
    z_dm <- NULL
    int_z_dm_keep_g <- NULL
    dat_reg[, z_dm := get(instrument) - mean(get(instrument))]
    dat_reg[, keep_g := ifelse(is.nan(t_val), FALSE, t_val >= qnorm(1-p_th))]
    dat_reg[, int_z_dm_keep_g := z_dm * keep_g]
    # f4 <- formula(paste0(yname, " ~ keep_g | ", treat, " ~ z_dm:keep_g"))
    f4 <- formula(paste0(yname, " ~ keep_g | ", treat, " ~ int_z_dm_keep_g"))
    reg <- fixest::feols(data = dat_reg, f4, vcov = "hetero")
    return(reg)
  }
}


#' Create groups from score prediction
#'
#' @import data.table
#' @import grf
#' @import gtools
#'
#' @param data A data.frame containing the necessary variables to run the model.
#' @param treat Name of the treatment variable
#' @param instrument Name of the instrument
#' @param controls Controls. Either as a formula, eg ~ x1 + x2, or a vector of strings, eg c("x1", "x2")
#' @param n_groups Number of groups to cut score into
#' @param pred_method Score prediction method
#' @param n_folds Number of folds to use for score prediction
#'
#' @return data.table
#' @export
#'
#' @examples
#' # Loading packages
#' library(fixest)
#'
#' # Generate a simulated dataset
#' data <- sim.data()
#'
#' data2 <- create.score.groups(data, treat = "d", instrument = "z", controls = ~x,
#'                              n_groups = 10)
#'
#' # Print summary of predicted compliance scores
#' summary(data2$pred_p)
#'
#' # Run a standard IV
#' reg1 <- feols(data = data2, y ~ 1 | d ~ z)
#'
#' # Run basic test-and-select method
#' reg2 <- run.late.rest(data = data2, yname = "y", treat = "d", instrument = "z", controls = ~g)
#'
#' # Run basic test-and-select method
#' reg3 <- run.late.rest(data = data2, yname = "y", treat = "d", instrument = "z",
#'                                    controls = ~score_g)
#'
#' # Comparing both regressions
#' etable(reg1, reg2, reg3)

create.score.groups <- function(data, treat, instrument, controls, n_groups,
                                pred_method = "Causal_Forest",
                                n_folds = 2){
  # Convert data to data.table
  dat <- as.data.table(data)

  # Load helper function ----------------------------------------------------
  rand_n_list_group <- function(n, n_groups){
    rem <- n %% n_groups
    mult <- n %/% n_groups
    tmp <- sample(rep(seq(1, n_groups), mult), mult*n_groups, replace = FALSE)
    return(c(tmp, sample(seq(1, n_groups), rem)))
  }

  # Parsing -----------------------------------------------------------------

  if(inherits(controls, "character")){
    model_str <- formula(paste("~-1+", paste(controls, collapse = "+")))
    cov_list <- controls
  }
  if(inherits(controls, "formula")){
    model_str <- update(controls, ~ -1 + . )
    cov_list <- all.vars(controls[[2]])
  }

  # Checks ------------------------------------------------------------------

  # Remove covariates with NAs
  any_na <- dat[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = cov_list]

  if(any(any_na > 0)){
    warning("Some of the control variables had missing values. Corresponding rows were removed.")
    dat <- dat[complete.cases(dat[, cov_list, with = FALSE])]
  }
  rm(any_na)

  # Checking group size
  if(nrow(dat)/n_groups <= 10){
    warning("The chosen number of groups will yield groups with less than 10 observations.")
  }

  # Code --------------------------------------------------------------------
  split_x = pred_p = score_g = NULL

  # Data split
  dat[, split_x := rand_n_list_group(.N, n_folds)]

  # Cross-fit heterogenous compliance score prediction
  dat[, pred_p := NA_real_]
  for(i in 1:n_folds){
    if(pred_method == "Causal_Forest"){
      pred_model <- causal_forest(
        model.matrix(model_str, dat[split_x != i]),
        matrix(dat[split_x != i][[treat]]),
        matrix(dat[split_x != i][[instrument]])
      )
      dat[split_x == i, pred_p := predict(pred_model,
                                           model.matrix(model_str, dat[split_x == i]))]
    } else {
      stop("Prediction method not recognized.")
    }
  }

  # Cut predicted compliance score into n_groups quantiles
  tryCatch(
    {
      dat[, score_g := quantcut(pred_p, q = n_groups, labels = FALSE)]
      tmp <- dat[, list(n = .N), by = score_g]
      tmp <- max(tmp$n) - min(tmp$n)
      if(tmp > 5){
        warning("Added noise of the order of 2^-30 to break score ties. Consider using less groups.")
        dat[, score_g := quantcut(pred_p + rnorm(.N, 0, 2^-30), q = n_groups, labels = FALSE)]
      }
      rm(tmp)
    },
    error = function(cond){
      warning("Added noise of the order of 2^-30 to break score ties. Consider using less groups.")
      dat[, score_g := quantcut(pred_p + rnorm(.N, 0, 2^-30), q = n_groups, labels = FALSE)]
    }
  )


  # Delete irrelevant data
  dat[, split_x := NULL]

  # Return data set
  return(dat)
}


#' D-LATE estimation with cross-fitted score estimation
#'
#' @import fixest
#' @import data.table
#' @import lmtest
#' @import sandwich
#' @import grf
#' @importFrom stats lm quantile qnorm
#'
#' @param data A data.frame containing the necessary variables to run the model.
#' @param yname Name of the dependent variable (character string)
#' @param treat Name of the treatment variable (character string)
#' @param instrument Name of the instrument (character string)
#' @param controls Controls. Either as a formula, eg ~ x1 + x2, or a vector of strings, eg c("x1", "x2")
#' @param n_groups Number of groups to split the predicted compliance scores into (default = 10)
#' @param pred_method Method for predicting compliance scores. Currently only "Causal_Forest" is supported
#' @param p_th Group selection threshold for the t-test (default = 0.05)
#' @param max_tries Maximum number of attempts to create a valid split with sufficient instrument variation (default = 3)
#' @param verbose Logical indicating whether to show warnings from t-statistic computations (default = FALSE)
#' @param weighted Logical indicating whether to use weighted IV rather than interacted IV (default = FALSE)
#'
#' @details
#' This function implements the D-LATE estimator with the following steps:
#' 1. Splits the data into 3 folds
#' 2. For each combination of folds:
#'    - Uses one fold to predict compliance scores and determine group thresholds
#'    - Uses the second fold to compute first-stage t-statistics for each group
#'    - Assigns these t-statistics to observations in the third fold
#' 3. Performs final IV estimation using group selection based on t-statistics
#'
#' @return A fixest object containing the IV regression results
#' @export
#'
#' @examples
#' # Loading required packages
#' library(fixest)
#' library(grf)
#'
#' # Generate a simulated dataset
#' data <- sim.data()
#'
#' # Run standard IV
#' reg1 <- feols(data = data, y ~ 1 | d ~ z)
#'
#' # Run D-LATE
#' reg2 <- dlate_method(data = data,
#'                     yname = "y",
#'                     treat = "d",
#'                     instrument = "z",
#'                     controls = ~x)
#'
#' # Compare results
#' etable(reg1, reg2)

dlate_method <- function(data, yname, treat, instrument, controls,
                             n_groups = 10, pred_method = "Causal_Forest",
                             p_th = 0.05, max_tries = 3, verbose = FALSE,
                         weighted = FALSE) {
  # Input checks
  if (!inherits(data, "data.frame")) {
    stop("'data' must be a data.frame or data.table")
  }
  if (!is.character(yname) || length(yname) != 1) {
    stop("'yname' must be a single character string")
  }
  if (!is.character(treat) || length(treat) != 1) {
    stop("'treat' must be a single character string")
  }
  if (!all(unique(data[[treat]]) %in% c(0, 1))) {
    stop("Treatment variable must be binary (0/1)")
  }
  if (!is.character(instrument) || length(instrument) != 1) {
    stop("'instrument' must be a single character string")
  }
  if (!all(unique(data[[instrument]]) %in% c(0, 1))) {
    stop("Instrument variable must be binary (0/1)")
  }
  if (!(inherits(controls, "formula") || inherits(controls, "character"))) {
    stop("'controls' must be either a formula or character vector")
  }

  # Convert data to data.table
  dat <- as.data.table(data)

  # Parse controls
  if(inherits(controls, "character")) {
    model_str <- formula(paste("~-1+", paste(controls, collapse = "+")))
    cov_list <- controls
  }
  if(inherits(controls, "formula")) {
    model_str <- update(controls, ~ -1 + . )
    cov_list <- all.vars(controls[[2]])
  }

  # Check for and handle NAs
  any_na <- dat[, lapply(.SD, function(x) sum(is.na(x))),
                .SDcols = c(yname, treat, instrument, cov_list)]

  if(any(unlist(any_na) > 0)) {
    warning("Some variables had missing values. Corresponding rows were removed.")
    dat <- dat[complete.cases(dat[, c(yname, treat, instrument, cov_list), with = FALSE])]
  }

  # Check group size
  if(nrow(dat)/n_groups <= 10) {
    warning("The chosen number of groups will yield groups with less than 10 observations.")
  }

  # Create helper function for random splits
  rand_n_list_group <- function(n, n_groups) {
    rem <- n %% n_groups
    mult <- n %/% n_groups
    tmp <- sample(rep(seq(1, n_groups), mult), mult*n_groups, replace = FALSE)
    return(c(tmp, sample(seq(1, n_groups), rem)))
  }

  # Function to get quantile thresholds
  get_quantile_thresholds <- function(x, n_groups) {
    probs <- seq(0, 1, length.out = n_groups + 1)
    unique(quantile(x, probs = probs, type = 1))
  }

  # Variables for R check
  split_x = pred_p = score_g = t_val = z_dm = keep_g = NULL

  for(try_num in 1:max_tries) {

    # Initialize prediction column
    dat[, split_x := NA_integer_]
    dat[, pred_p := NA_real_]
    dat[, score_g := NA_integer_]
    dat[, t_val := NA_real_]

    # Create three folds
    dat[, split_x := rand_n_list_group(.N, 3)]

    # For each combination of folds
    fold_combos <- list(
      list(pred = 1, fs = 2, final = 3),
      list(pred = 2, fs = 3, final = 1),
      list(pred = 3, fs = 1, final = 2)
    )

    tryCatch({
      for(combo in fold_combos) {
        # Predict first stage in prediction fold
        if(pred_method == "Causal_Forest") {
          pred_model <- causal_forest(
            model.matrix(model_str, dat[split_x == combo$pred]),
            matrix(dat[split_x == combo$pred][[treat]]),
            matrix(dat[split_x == combo$pred][[instrument]])
          )

          # Get thresholds from prediction fold
          pred_fold_preds <- predict(pred_model,
                                     model.matrix(model_str, dat[split_x == combo$pred]))$predictions
          thresholds <- get_quantile_thresholds(pred_fold_preds, n_groups)

          # Check if we have fewer unique thresholds than requested groups
          n_actual_groups <- length(thresholds) - 1
          if(n_actual_groups < n_groups) {
            warning("Due to ties in predictions, only ", n_actual_groups,
                    " distinct groups can be created instead of the requested ", n_groups, " groups")
          }

          # Apply to other folds
          for(fold in c(combo$fs, combo$final)) {
            fold_preds <- predict(pred_model,
                                  model.matrix(model_str, dat[split_x == fold]))
            dat[split_x == fold, pred_p := fold_preds]
            dat[split_x == fold, score_g := cut(pred_p,
                                                breaks = thresholds,
                                                labels = FALSE,
                                                include.lowest = TRUE)]
          }
        } else {
          stop("Prediction method not recognized.")
        }

        # Compute first-stage t-values in FS fold and save for final fold
        # f1 <- formula(paste0(treat, " ~ ", instrument))
        if (!verbose) {
          t_vals <- suppressWarnings(
            dat[split_x == combo$fs,
                        as.list(lmtest::coeftest(
                          lm(.SD[[1]] ~ .SD[[2]], data = .SD),
                          vcov = sandwich::vcovHC,
                          type = "HC3")[2, c(1, 3)]),
                        .SDcols = c(treat, instrument),
                        by = score_g]
          )
        } else {
          # Allow warnings to show
          t_vals <- dat[split_x == combo$fs,
                        as.list(lmtest::coeftest(
                          lm(.SD[[1]] ~ .SD[[2]], data = .SD),
                          vcov = sandwich::vcovHC,
                          type = "HC3")[2, c(1, 3)]),
                        .SDcols = c(treat, instrument),
                        by = score_g]
        }
        names(t_vals)[2:3] <- c("pc", "t_val")

        # Merge t-values to final fold
        dat[split_x == combo$final, t_val := t_vals[.SD, on = "score_g"]$t_val]
      }

      # If we get here, everything worked
      if(try_num > 1) {
        warning("Had to retry split creation ", try_num - 1, " time(s) due to constant instruments within groups")
      }
      break  # Exit the retry loop

    }, error = function(e) {
      if(try_num == max_tries) {
        stop("Failed to create valid split after ", max_tries, " attempts. Check instrument variation in your data")
      }
      # Clear any partial results
      dat[, c("pred_p", "score_g", "t_val") := list(NA_real_, NA_integer_, NA_real_)]
    })
  }

  # Prepare final regression
  dat[, z_dm := .SD[[1]] - mean(.SD[[1]]), .SDcols = instrument]
  # dat[, keep_g := ifelse(is.nan(t_val), FALSE, t_val >= qnorm(1-p_th))]
  dat[, keep_g := as.numeric(ifelse(is.nan(t_val), FALSE, t_val >= qnorm(1-p_th)))]
  int_zdm_keep_g <- NULL
  dat[, int_zdm_keep_g := z_dm * keep_g]

  # Run final regression
  if(weighted == FALSE){
    # f4 <- formula(paste0(yname, " ~ keep_g | ", treat, " ~ z_dm:keep_g"))
    f4 <- formula(paste0(yname, " ~ keep_g | ", treat, " ~ int_zdm_keep_g"))
    reg <- fixest::feols(data = dat, f4, vcov = "hetero")
  }
  if(weighted == TRUE){
    f4 <- formula(paste0(yname, " ~ 1 | ", treat, " ~ z"))
    reg <- fixest::feols(data = dat, f4, vcov = "hetero", weights = ~keep_g)
  }

  return(reg)
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
