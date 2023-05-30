run_bart <- function(file, file_indep, file_ref) {
  on.exit({rm(list=ls()); gc()})
  dat <- read.csv(file)
  dat_indep <- read.csv(file_indep)
  dat_ref <- read.csv(file_ref)
  colnames(dat_indep) <- paste0(colnames(dat_indep), "_indep")
  colnames(dat_ref) <- paste0(colnames(dat_ref), "_ref")
  dat <- dat %>%
    bind_cols(dat_indep, dat_ref) %>%
    mutate(obs = 1:n())
  x_mat <- as.matrix(select(dat, matches("X[[:digit:]]*$")))
  post <- surv.bart(
    x.train = x_mat,
    times = dat$Y_1,
    delta = dat$Y_2,
    K = 100,
    nskip = 1e4,
    ndpost = 1e4,
    keepevery = 1,
    printevery = 1e10
  )
  times <- c(.11, .36, 1.24, 4.08, 10.36)
  x_pred <- expand_grid(t = times, dat)
  # pred <- predict(post, newdata = as.matrix(select(x_pred, t, x1, x2)))
  pred_indep <- predict(
    post,
    newdata = as.matrix(select(x_pred, t, matches("X[[:digit:]]*_indep")))
  )
  pred_ref <- predict(
    post,
    newdata = as.matrix(select(x_pred, t, matches("X[[:digit:]]*_ref")))
  )

  beta <- 2
  betas <- c(-1, 1, 2, 0, 0, -2, -1, 1, 1.5, -1.5)
  x_ref <- as.matrix(select(x_pred, matches("X[[:digit:]]*_ref")))
  x_indep <- as.matrix(select(x_pred, matches("X[[:digit:]]*_indep")))
  hazard_ref <- 1 / beta * exp(c(x_ref %*% betas))
  hazard_indep <- 1 / beta * exp(c(x_indep %*% betas))
  
  stats <- x_pred %>%
    mutate(
      mean_indep = pred_indep$surv.test.mean,
      mean_ref = pred_ref$surv.test.mean,
      lb_ref = apply(pred_ref$surv.test, 2, quantile, prob = .025),
      ub_ref = apply(pred_ref$surv.test, 2, quantile, prob = .975),
      surv_true_ref = exp(-!!hazard_ref * t),
      surv_true_indep = exp(-!!hazard_indep * t),
      bias = mean_ref - surv_true_ref,
      coverage = lb_ref < surv_true_ref & surv_true_ref < ub_ref,
      brier = (as.numeric(Y_true_indep > t) - mean_indep) ^ 2,
      brier_true = (as.numeric(Y_true_indep > t) - surv_true_indep) ^ 2
    )
  brier <- stats %>%
    group_by(t) %>%
    summarize(
      brier = mean(brier),
      brier_true = mean(brier_true),
      .groups = "drop"
    )
  
  select(stats, t, obs, bias, coverage) %>%
    left_join(brier, by = "t")
}

get_ocs <- function(x) {
  x %>%
    group_by(t) %>%
    summarize(
      bias = mean(bias),
      rmse = sqrt(mean(sum(bias ^ 2))),
      coverage = mean(coverage),
      brier = mean(brier)
    )
}
