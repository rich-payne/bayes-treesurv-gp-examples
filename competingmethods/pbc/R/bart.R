fit_bart <- function(data, k_fold = NULL) {
  data <- get_k_fold_data(data, k_fold)
  x_mat <- get_x_mat(data)
  post <- surv.bart(
    x.train = x_mat,
    times = data$time,
    delta = data$status,
    nskip = 1e3,
    ndpost = 1e4,
    keepevery = 1,
    printevery = 1e10
  )
}

get_x_mat <- function(data) {
  x_mat <- model.matrix(
    ~ -1 + trt + age + sex + edema + bili + albumin + platelet,
    data = data
  )
}
