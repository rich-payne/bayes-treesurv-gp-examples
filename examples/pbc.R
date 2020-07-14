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
