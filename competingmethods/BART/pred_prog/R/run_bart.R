run_bart <- function(file, file_indep, file_ref, sim) {
  on.exit({rm(list=ls()); gc()})
  dat <- read.csv(file)
  dat_indep_orig <- dat_indep <- read.csv(file_indep)
  dat_ref_orig <- dat_ref <- read.csv(file_ref)
  colnames(dat_indep) <- paste0(colnames(dat_indep), "_indep")
  colnames(dat_ref) <- paste0(colnames(dat_ref), "_ref")
  dat <- dat %>%
    bind_cols(dat_indep, dat_ref) %>%
    mutate(obs = 1:n())
  x_mat <- as.matrix(select(dat, biomarker, trt))
  post <- surv.bart(
    x.train = x_mat,
    times = dat$Y_1,
    delta = dat$Y_2,
    K = 100,
    nskip = 1e3,
    ndpost = 1e4,
    keepevery = 1,
    printevery = 1e10
  )
  times <- c(.55, .81, 1.13, 1.84, 2.50)
  x_pred <- expand_grid(dat, t = times)
  # pred <- predict(post, newdata = as.matrix(select(x_pred, t, x1, x2)))
  # pred_indep <- predict(
  #   post,
  #   newdata = as.matrix(select(x_pred, t, biomarker_indep, trt_indep))
  # )
  # pred_ref <- predict(
  #   post,
  #   newdata = as.matrix(select(x_pred, t, biomarker_ref, trt_ref))
  # )
  pred_indep <- pred_bart(post, dat_indep_orig, times)
  pred_ref <- pred_bart(post, dat_ref_orig, times)

  a_shape <- 1
  b_shape <- 5
  a_scale <- 1
  b_scale <- 2
  
  g <- get_cens_dist(dat)
  
  stats <- x_pred %>%
    mutate(
      mean_indep = pred_indep$surv,
      mean_ref = pred_ref$surv,
      lb_ref = pred_ref$lb,
      ub_ref = pred_ref$ub,
      w_shape_ref = a_shape + b_shape * biomarker_ref,
      w_scale_ref = a_scale + b_scale * trt_ref * biomarker_ref,
      w_shape_indep = a_shape + b_shape * biomarker_indep,
      w_scale_indep = a_scale + b_scale * trt_indep * biomarker_indep,
      surv_true_ref = pweibull(t, w_shape_ref, w_scale_ref, lower = FALSE),
      surv_true_indep = pweibull(t, w_shape_indep, w_scale_indep, lower = FALSE),
      bias = mean_ref - surv_true_ref,
      coverage = lb_ref < surv_true_ref & surv_true_ref < ub_ref,
      # calculate brier
      val1 = mean_indep ^ 2,
      val2 = (1 - mean_indep) ^ 2,
      w1 = safe_mult(1 / g(Y_1_indep), i(Y_1_indep <= t) * i(Y_2_indep == 1)),
      w2 = safe_mult(1 / g(t), i(Y_1_indep > t)),
      brier_cens = w1 * val1 + w2 * val2,
      brier = (as.numeric(Y_true_indep > t) - mean_indep) ^ 2,
      brier_true = (as.numeric(Y_true_indep > t) - surv_true_indep) ^ 2
    )
  stats <- select(stats, t, obs, bias, coverage, starts_with("brier")) %>%
    mutate(sim = !!sim)
  return(stats)
}

get_cens_dist <- function(data) {
  fit <- survfit(Surv(Y_1, abs(1 - Y_2)) ~ 1, data = data)
  # censoring distribution
  g <- approxfun(
    fit$time,
    fit$surv,
    method = "constant",
    f = 0,
    yleft = 1,
    rule = 2
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

pred_bart <- function(post, data, times) {
  lb <- ub <- surv <- matrix(NA, length(times), nrow(data))
  for (i in 1:nrow(data)) {
    xmat <- select(data, biomarker, trt) %>%
      slice(i) %>%
      expand_grid(t = !!post$times) %>%
      select(t, everything()) %>%
      as.matrix()
    pred <- predict(post, newdata = xmat)
    lower_bnd <- approx(
      post$times,
      apply(pred$surv.test, 2, quantile, p = .025), times
    )$y
    upper_bnd <- approx(
      post$times,
      apply(pred$surv.test, 2, quantile, p = .975), times
    )$y
    surv[, i] <- approx(post$times, pred$surv.test.mean, times)$y
    lb[, i] <- lower_bnd
    ub[, i] <- upper_bnd
  }
  data.frame(
    surv = c(surv),
    lb = c(lb),
    ub = c(ub)
  )
}

add_brier_miss <- function(barts, missing_data) {
  u_times <- sort(unique(barts$t))
  miss_long <- missing_data %>%
    group_by(sim) %>%
    mutate(obs = 1:n()) %>%
    ungroup() %>%
    pivot_longer(
      -c("obs", "sim"),
      names_to = "t",
      values_to = "missing",
      names_pattern = "time([[:digit:]])",
      names_transform = list(t = as.numeric)
    ) %>%
    mutate(t = u_times[t])
  brier_miss <- barts %>%
    left_join(miss_long, by = c("sim", "t", "obs")) %>%
    filter(missing == 0) %>%
    group_by(t, sim) %>%
    summarize(
      brier_cens_miss = mean(brier_cens),
      .groups = "drop"
    )
  brier <- barts %>%
    group_by(t, sim) %>%
    summarize(
      brier = mean(brier),
      brier_cens = mean(brier_cens),
      .groups = "drop"
    ) %>%
    left_join(brier_miss, by = c("t", "sim"))
}
