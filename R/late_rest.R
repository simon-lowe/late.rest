#' Main function
#'
#' @import fixest
#' @import data.table
#' @import lmtest
#' @import sandwich
#' @import readr
#' @importFrom stats lm
#'
#' @param data data.frame
#' @param yname dependent variable name
#' @param treat bla
#' @param instrument blu
#' @param controls shu
#' @param n_folds sha
#' @param p_th Group selection threshold
#'
#' @return tada
#' @export
#'
#' @examples
#' # You can do the following hihi

run.late.rest <- function(data, yname, treat, instrument, controls, n_folds = 2, p_th = 0.05){

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

  # # Create data
  # data <- dat_reg
  # setnames(dat_reg, c(yname, treat, instrument), c("y", "d", "z"))

  # Create the group variable
  data[, g := .GRP, by = mget(c_vars)]

  # Create the split
  data[, split_x := rand_n_list_group(.N, n_folds), by = g]

  # return(c(min(data[, .(n = .N), by = g]$n), data[, .(n = .N), by = split_x]$n))

  # Compute the cross-fitted compliance rates

  f1 <- formula(paste0(treat, " ~ ", instrument))

  tmp <- data[, lapply(1:n_folds, function(x) {
    # bla <- lmtest::coeftest(lm(d ~ z, subset = split_x != x), vcov = sandwich::vcovHC, type = "HC3")
    bla <- lmtest::coeftest(lm(data = .SD, formula = f1, subset = split_x != x), vcov = sandwich::vcovHC, type = "HC3")
    list(bla[2,1], bla[2, 4])
  }), keyby = g]
  tmp[, res := rep(c("pc", "p_val"), .N/2)]
  # tmp <- data.table::melt(tmp, measure = data.table:::patterns("V\\d+"), value.name = "val", variable.name = "split_x", variable.factor = F)
  tmp <- data.table::melt(tmp, measure = grep("^V\\d+$", colnames(tmp), value = TRUE), value.name = "val", variable.name = "split_x", variable.factor = F)
  tmp[, split_x := readr::parse_number(split_x)]
  tmp[, val := unlist(val)]
  tmp <- data.table::dcast(tmp, g + split_x ~ res, value.var = "val")

  dat_reg <- data.table::merge.data.table(data, tmp, by = c("g", "split_x"))

  f2 <- formula(paste0(yname, " ~ 1 | ", treat, " ~ ", instrument))

  reg <- fixest::feols(data = dat_reg[p_val <= p_th], f2)

  return(reg)
}
