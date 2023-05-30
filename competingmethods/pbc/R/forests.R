fit_forest <- function(data, k_fold = NULL) {
  data <- get_k_fold_data(data, k_fold)
  fit <- rfsrc(
    Surv(time, status) ~ trt + age + sex + edema + bili + albumin + platelet,
    data = data,
    ntree = 1000,
    block.size = 1
  )
}

get_k_fold_data <- function(data, k_fold) {
  if (!is.null(k_fold)) {
    data <- filter(data, k_fold != !!k_fold)
  }
  return(data)
}