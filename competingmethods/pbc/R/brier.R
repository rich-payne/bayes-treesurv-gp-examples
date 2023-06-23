get_brier_forest <- function(forest_fit, times, data, fold, missing_data) {
  data_fold <- filter(data, k_fold == !!fold)
  missing_data_fold <- missing_data %>% filter(fold == !!fold) %>%
    select(-fold) %>%
    as.matrix()
  g <- get_cens_dist(data)
  pred <- pred_forest(forest_fit, newdata = data_fold, times = times)
  brier <- purrr::imap_dfr(
    times,
    get_brier,
    pred = pred,
    y = data_fold$time,
    delta = data_fold$status,
    g = g,
    missing_data_fold = missing_data_fold
  ) %>%
    mutate(fold = !!fold, model = "forest")
}

get_brier_bart <- function(bart_fit, times, data, fold, missing_data) {
  data_fold <- filter(data, k_fold == !!fold)
  missing_data_fold <- missing_data %>% filter(fold == !!fold) %>%
    select(-fold) %>%
    as.matrix()
  g <- get_cens_dist(data)
  x_mat <- get_x_mat(data)
  x_mat <- x_mat[data$k_fold == fold, ]
  pred <- pred_bart(bart_fit, newdata = x_mat, times = times)
  brier <- purrr::imap_dfr(
    times,
    get_brier,
    pred = pred,
    y = data_fold$time,
    delta = data_fold$status,
    g = g,
    missing_data_fold = missing_data_fold
  ) %>%
    mutate(fold = !!fold, model = "BART")
}

get_cens_dist <- function(data) {
  fit <- survfit(Surv(time, abs(1 - status)) ~ 1, data = data)
  # censoring distribution
  g <- approxfun(
    fit$time,
    fit$surv,
    method = "constant",
    f = 0
  )
  return(g)
}

i <- function(x) as.numeric(x)

safe_mult <- function(x, y) {
  ind <- which(x == Inf & y == 0)
  out <- x * y
  out[ind] <- 0
  return(out)
}

pred_forest <- function(forest_fit, newdata, times) {
  pred <- predict(forest_fit, newdata = newdata)
  pred <- apply(
    pred$survival,
    1,
    function(surv) approx(pred$time.interest, surv, times)$y
  )
}

pred_bart <- function(bart_fit, newdata, times) {
  surv <- matrix(NA, length(times), nrow(newdata))
  newdata <- as_tibble(newdata)
  for (i in 1:nrow(newdata)) {
    xmat <- newdata %>%
      slice(i) %>%
      expand_grid(t = !!bart_fit$times) %>%
      select(t, everything()) %>%
      as.matrix()
    pred <- predict(bart_fit, newdata = xmat)
    surv[, i] <- approx(bart_fit$times, pred$surv.test.mean, times)$y
  }
  return(surv)
}

get_brier <- function(time, ind, pred, y, delta, g, missing_data_fold = NULL) {
  # survival function at specific time for all subjects,
  #   i.e. times x subjects
  if (is.null(missing_data_fold)) {
    missing_data_fold <- rep(0, ncol(pred))
  }
  ind_miss <- missing_data_fold[, ind] == 1
  
  surv <- pred[ind, ]
  val1 <- surv ^ 2
  val2 <- (1 - surv) ^ 2
  w1 <- safe_mult(1 / g(y), i(y <= time) * i(delta == 1))
  w2 <- safe_mult(1 / g(time), i(y > time))
  brier <- w1 * val1 + w2 * val2
  data.frame(
    time = time,
    brier_miss = mean(brier[!ind_miss]),
    brier = mean(brier)
  )
}
