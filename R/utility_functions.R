rand_n_list_group <- function(n, n_groups){
  rem <- n %% n_groups
  mult <- n %/% n_groups
  tmp <- sample(rep(seq(1, n_groups), mult), mult*n_groups, replace = FALSE)
  return(c(tmp, sample(seq(1, n_groups), rem)))
}
