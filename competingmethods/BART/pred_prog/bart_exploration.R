library(targets)
library(dplyr)
library(tidyr)
library(BART)
bart <- readRDS("bart_204.rds")
tar_load(c(file_204, file_104, data_file_ref))
data <- read.csv(file_204)
data_indep <- read.csv(file_104)
data_ref <- read.csv(data_file_ref)

times <- bart$fit$times
xmat <- select(data, biomarker, trt) %>%
  slice_head(n = 1) %>%
  expand_grid(t = !!times) %>%
  select(t, everything()) %>%
  as.matrix()
pred <- predict(bart$fit, newdata = xmat)
plot(pred$times, pred$surv.test.mean, type = "l")

a_shape <- 1
b_shape <- 5
a_scale <- 1
b_scale <- 2

w_shape <- a_shape + b_shape * xmat[, "biomarker"]
w_scale <- a_scale + b_scale * xmat[, "trt"] * xmat[, "biomarker"]
surv_true <- pweibull(xmat[, "t"], w_shape, w_scale, lower = FALSE)
lines(xmat[, "t"], surv_true, col = "red")

# # try a hold out dataset (independent)
# xmat_indep <- select(data_indep, biomarker, trt) %>%
#   slice_head(n = 1) %>%
#   expand_grid(t = !!times) %>%
#   # mutate(t = 1) %>%
#   select(t, everything()) %>%
#   as.matrix()
# pred_indep <- predict(bart$fit, newdata = xmat_indep)
# plot(pred_indep$times, pred_indep$surv.test.mean, type = "l")
# 
# # subset of times?
# # try a hold out dataset (independent)
# xmat_indep_sub <- select(data_indep, biomarker, trt) %>%
#   slice_head(n = 1) %>%
#   expand_grid(t = 1:2) %>%
#   # mutate(t = 1) %>%
#   select(t, everything()) %>%
#   as.matrix()
# pred_indep_sub <- predict(bart$fit, newdata = xmat_indep_sub)
# points(pred_indep_sub$times, pred_indep_sub$surv.test.mean)
# 
# # with pre.bart?
# pre <- surv.pre.bart(
#   times = data$Y_1,
#   delta = data$Y_2,
#   x.train = as.matrix(select(data, biomarker, trt)),
#   x.test = as.matrix(select(data, biomarker, trt) %>% slice_head(n = 1))
# )
# 
# w_shape_indep <- a_shape + b_shape * xmat_indep[, "biomarker"]
# w_scale_indep <- a_scale + b_scale * xmat_indep[, "trt"] * xmat_indep[, "biomarker"]
# surv_true_indep <- pweibull(xmat[, "t"], w_shape_indep, w_scale_indep, lower = FALSE)
# lines(xmat_indep[, "t"], surv_true_indep, col = "red")

# calculate a brier score

qtimes <- quantile(data$Y_1, c(0.1, 0.25, 0.5, 0.75, 0.90))
briers <- matrix(NA, nrow(data_indep), length(qtimes))
for (i in 1:nrow(data_indep)) {
  print(i)
  xmat_indep <- select(data_indep, biomarker, trt) %>%
    slice(i) %>%
    expand_grid(t = !!times) %>%
    select(t, everything()) %>%
    as.matrix()
  pred <- predict(bart$fit, newdata = xmat_indep)
  surv_true_indep <- approx(times, pred$surv.test.mean, qtimes)$y
  briers[i, ] <- (as.numeric(data_indep$Y_true[i] > qtimes) - surv_true_indep) ^ 2
}
plot(qtimes, apply(briers, 2, mean))
