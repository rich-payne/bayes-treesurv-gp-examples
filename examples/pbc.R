library(survival)
library(dplyr)
out <- pbc %>%
  dplyr::filter(trt %in% c(1, 2), complete.cases(.)) %>%
  dplyr::mutate(
    status = case_when(
      status == 1 ~ 0L,
      status == 2 ~ 1L,
      TRUE ~ status
    ),
    trt = case_when(
      trt == 2 ~ 0L,
      TRUE ~ trt
    )
  )
write.csv(out, "pbc.csv", row.names = FALSE)

group_by(out, trt) %>%
  summarize(
    mean_time = mean(time),
    sd_time = sd(time)
  )

# split for k-fold cross validation
set.seed(883311)
k <- 10
out_k <- out %>%
  slice_sample(prop = 1) %>%
  mutate(k_fold = cut(1:n(), breaks = !!k, label = FALSE))
write.csv(out_k, "pbc_kfold.csv", row.names = FALSE)


