#' Generate an group index list of length n with n_groups groups
#'
#' @param n Length of the list
#' @param n_groups Number of groups
#'
#' @return An integer vector
#' @export
#'
#' @examples
#' # Group index list with 10 observations and 3 groups
#' rand_n_list_group(10, 3)

rand_n_list_group <- function(n, n_groups){
  rem <- n %% n_groups
  mult <- n %/% n_groups
  tmp <- sample(rep(seq(1, n_groups), mult), mult*n_groups, replace = FALSE)
  return(c(tmp, sample(seq(1, n_groups), rem)))
}
