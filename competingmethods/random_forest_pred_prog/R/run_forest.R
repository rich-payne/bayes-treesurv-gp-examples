run_forest <- function(file, file_indep, file_ref, sim) {
  on.exit({rm(list=ls()); gc()})
  dat <- read.csv(file)
  dat_indep <- read.csv(file_indep)
  dat_ref <- read.csv(file_ref)
  colnames(dat_indep) <- paste0(colnames(dat_indep), "_indep")
  colnames(dat_ref) <- paste0(colnames(dat_ref), "_ref")
  dat <- dat %>%
    bind_cols(dat_indep, dat_ref) %>%
    mutate(obs = 1:n())
  fit <- rfsrc(
    Surv(Y_1, Y_2) ~ biomarker + trt,
    data = dat,
    ntree = 1000,
    block.size = 1
  )
  times <- c(.55, .81, 1.13, 1.84, 2.50)
  # get close enough.
  approx_times_index <- vapply(
    times,
    function(time, model_times) {
      which.min(abs(time - model_times))
    },
    numeric(1),
    model_times = fit$time.interest
  )
  approx_times <- fit$time.interest[approx_times_index]
  
  pred_ref <- get_forest_pred(
    fit,
    select(dat, biomarker = biomarker_ref, trt = trt_ref),
    approx_times_index
  )
  pred_indep <- get_forest_pred(
    fit,
    select(dat, biomarker = biomarker_indep, trt = trt_indep),
    approx_times_index
  )
  
  sumry_ref <- pred_ref %>%
    group_by(time, data_index) %>%
    summarize(
      mean_ref = mean(survival),
      lb_ref = quantile(survival, prob = .025),
      ub_ref = quantile(survival, prob = .975),
      .groups = "drop"
    )
  sumry_indep <- pred_indep %>%
    group_by(time, data_index) %>%
    summarize(
      mean_indep = mean(survival),
      lb_indep = quantile(survival, prob = .025),
      ub_indp = quantile(survival, prob = .975),
      .groups = "drop"
    )
  sumry_all <- left_join(sumry_indep, sumry_ref, by = c("time", "data_index"))
  
  x_pred <- expand_grid(
    data.frame(t = approx_times, true_time = times),
    mutate(dat, data_index = 1:n())
  )
  
  a_shape <- 1
  b_shape <- 5
  a_scale <- 1
  b_scale <- 2
  
  g <- get_cens_dist(dat)
  
  stats <- x_pred %>%
    left_join(sumry_all, by = c(t = "time", "data_index")) %>%
    mutate(
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
  brier <- select(stats, time = t, data_index, starts_with("brier")) %>%
    mutate(sim = sim)
  brier_sumry <- stats %>%
    group_by(true_time) %>%
    summarize(
      brier = mean(brier),
      brier_cens = mean(brier_cens),
      brier_true = mean(brier_true),
      .groups = "drop"
    ) %>%
    mutate(sim = sim)
  stats <- select(stats, t = true_time, obs, bias, coverage,starts_with("brier")) %>%
    mutate(sim = sim)
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

get_forest_pred <- function(fit, data, approx_times_index) {
  purrr::map_dfr(
    seq_len(fit$ntree),
    function(tree_num, data, fit, cols) {
      pred <- predict(
        fit,
        newdata = data,
        get.tree = tree_num
      )
      out <- data.frame(
        tree = tree_num,
        data_index = rep(seq_len(nrow(data)), length(cols)),
        time = rep(pred$time.interest[cols], each = nrow(data)),
        survival = c(pred$survival[, cols])
      ) 
    },
    data = data,
    fit = fit,
    cols = approx_times_index
  ) 
}

add_brier_miss <- function(forests, missing_data) {
  u_times <- sort(unique(forests$t))
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
  brier_miss <- forests %>%
    left_join(miss_long, by = c("sim", "t", "obs")) %>%
    filter(missing == 0) %>%
    group_by(t, sim) %>%
    summarize(
      brier_cens_miss = mean(brier_cens),
      .groups = "drop"
    )
  brier <- forests %>%
    group_by(t, sim) %>%
    summarize(
      brier = mean(brier),
      brier_cens = mean(brier_cens),
      .groups = "drop"
    ) %>%
    left_join(brier_miss, by = c("t", "sim"))
}
