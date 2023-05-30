run_bart <- function(file, file_indep, file_ref) {
  on.exit({rm(list=ls()); gc()})
  dat <- read.csv(file)
  dat_indep <- read.csv(file_indep)
  dat_ref <- read.csv(file_ref)
  colnames(dat_indep) <- paste0(colnames(dat_indep), "_indep")
  colnames(dat_ref) <- paste0(colnames(dat_ref), "_ref")
  dat <- dat %>%
    bind_cols(dat_indep, dat_ref) %>%
    mutate(
      x2 = as.integer(factor(x2)),
      x2_indep = as.integer(factor(x2_indep)),
      x2_ref = as.integer(factor(x2_ref)),
      obs = 1:n()
    )
  x_mat <- as.matrix(select(dat, x1, x2))
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
  times <- c(.09, .62, 2.04, 4.6, 6.46)
  x_pred <- expand_grid(t = times, dat)
  # pred <- predict(post, newdata = as.matrix(select(x_pred, t, x1, x2)))
  pred_indep <- predict(
    post,
    newdata = as.matrix(select(x_pred, t, x1 = x1_indep, x2 = x2_indep))
  )
  pred_ref <- predict(
    post,
    newdata = as.matrix(select(x_pred, t, x1 = x1_ref, x2 = x2_ref))
  )

  stats <- x_pred %>%
    # group_by(t) %>%
    # mutate(subject = 1:n()) %>%
    # ungroup() %>%
    mutate(
      mean_indep = pred_indep$surv.test.mean,
      mean_ref = pred_ref$surv.test.mean,
      lb_ref = apply(pred_ref$surv.test, 2, quantile, prob = .025),
      ub_ref = apply(pred_ref$surv.test, 2, quantile, prob = .975),
      surv_true_ref = case_when(
        x2_ref <= 2 & x1_ref > 5 ~ pweibull(t, 2, 5, lower.tail = FALSE),
        x2_ref == 1 & x1_ref <= 5 ~ pweibull(t, 5, 1, lower.tail = FALSE),
        x2_ref == 2 & x1_ref <= 5 ~ pweibull(t, .9, 0.5, lower.tail = FALSE),
        x2_ref > 2 & x1_ref <= 3 ~ pweibull(t, 5, 5, lower.tail = FALSE),
        x2_ref > 2 & x1_ref > 3 & x1_ref <= 7 ~ pweibull(t, .5, .5, lower.tail = FALSE),
        x2_ref > 2 & x1_ref > 7 ~ custom_surv(t),
      ),
      surv_true_indep = case_when(
        x2_indep <= 2 & x1_indep > 5 ~ pweibull(t, 2, 5, lower.tail = FALSE),
        x2_indep == 1 & x1_indep <= 5 ~ pweibull(t, 5, 1, lower.tail = FALSE),
        x2_indep == 2 & x1_indep <= 5 ~ pweibull(t, .9, 0.5, lower.tail = FALSE),
        x2_indep > 2 & x1_indep <= 3 ~ pweibull(t, 5, 5, lower.tail = FALSE),
        x2_indep > 2 & x1_indep > 3 & x1_indep <= 7 ~ pweibull(t, .5, .5, lower.tail = FALSE),
        x2_indep > 2 & x1_indep > 7 ~ custom_surv(t),
      ),
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

custom_surv <- function(t) {
  times <- 0:11
  surv <- c(1, .9, .85, .5, .45, .44, .43, .1, .09, .05, .01, 0)
  ht <- -log(surv)
  theta1 <- 1
  surv1 <- exp(-ht * theta1)
  approxfun(times, surv)(t)
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
