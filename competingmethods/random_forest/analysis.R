library(randomForestSRC)
library(dplyr)
library(purrr)
library(ggplot2)

data <- read.csv(
  "pred_prog.csv",
  header = FALSE,
  col.names = c("time", "observed", "biomarker", "trt")
)

fit <- rfsrc(
  Surv(time, observed) ~ .,
  data = data,
  ntree = 1000,
  block.size = 1
)

# get survival predictions for each tree and calculate quantiles
pred <- purrr::map_dfr(
  seq_len(fit$ntree),
  function(tree_num) {
    pred <- predict(
      fit,
      newdata = data.frame(trt = 1, biomarker = 0.10),
      get.tree = tree_num
    )
    data.frame(
      tree = tree_num,
      time = pred$time.interest,
      survival = pred$survival[1, ]
    )
  }
) %>%
  group_by(time) %>%
  summarize(
    mean = mean(survival),
    lb = quantile(survival, prob = .025, names = FALSE),
    ub = quantile(survival, prob = .975, names = FALSE),
    .groups = "drop"
  )

ggplot(pred, aes(time, mean)) +
  geom_line() +
  geom_line(aes(y = lb), linetype = "dashed") +
  geom_line(aes(y = ub), linetype = "dashed")

