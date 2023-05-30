library(targets)

tar_option_set(
  packages = c("randomForestSRC", "dplyr")
)

list(
  tar_target(data_file, "../../examples/pred_prog.csv", format = "file"),
  tar_target(
    data,
    read.csv(
      data_file,
      header = FALSE,
      col.names = c("time", "observed", "biomarker", "trt")
    )
  ),
  tar_target(
    fit,
    rfsrc(Surv(time, observed) ~ ., data = data, ntree = 1000, block.size = 1)
  ),
  tar_target(
    pred,
    purrr::map_dfr(
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
  ),
  tar_target(
    random_forest,
    {
      filename <- "random_forest_pred.csv"
      write.csv(pred, filename)
      filename
    }
  )
)
