library(BART)
library(drake)
library(ggplot2)
library(dplyr)
library(tidyr)

get_x <- function(dat) {
  x <- as.matrix(dat[, 3:4])
  return(x)
}

# function to run analysis
run_bart <- function(dat) {
  x <- get_x(dat)
  post <- surv.bart(
    x.train = x,
    times = dat$time,
    delta = dat$observed,
    nskip = 1e3,
    ndpost = 1e3,
    keepevery = 2
  )
  return(post)
}

run_predict <- function(post, dat) {
  x <- get_x(dat)
  ind <- 436
  pre <- surv.pre.bart(
    times = dat$time,
    delta = dat$observed,
    x.train = x,
    x.test = x[ind, , drop = FALSE]
  )
  pred <- predict(post, newdata = pre$tx.test)
  return(pred)
}

get_plot_data <- function(post, pred) {
  qtls <- apply(pred$surv.test,2,quantile,probs=c(.025,.975))
  b <- pred$tx.test[1, 'biomarker']
  trt <- pred$tx.test[1, 'trt']
  true_surv <- 1 - pweibull(
    sort(unique(pred$tx.test[, 't'])),
    1 + 5 * b,
    1 + 2 * b * trt
  )
  dat_plot <- data.frame(
    time = post$times,
    mean = pred$surv.test.mean,
    lb = qtls[1, ],
    ub = qtls[2, ],
    true_surv = true_surv
  )
}

get_plot <- function(dat_plot) {
  p <- ggplot(dat_plot, aes(time)) +
    geom_line(aes(y = mean), linetype = 3) +
    geom_line(aes(y = lb), linetype = 2) +
    geom_line(aes(y = ub), linetype = 2) +
    geom_line(aes(y = true_surv)) +
    ylab("survival") +
    ggtitle("BART") +
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

plan <- drake_plan(
  dat = read.csv(
    file_in("../../examples/pred_prog.csv"),
    header = FALSE,
    col.names = c("time", "observed", "biomarker", "trt")
  ),
  post = run_bart(dat),
  pred = run_predict(post, dat),
  plot_data = get_plot_data(post, pred),
  plot_data_file = write.csv(plot_data, file_out("bart_pred.csv"), row.names = FALSE),
  bart_plot = get_plot(plot_data)
)

make(plan)
