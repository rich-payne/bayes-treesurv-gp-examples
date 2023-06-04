library(dplyr)
library(tidyr)
library(ggplot2)
library(fs)
library(purrr)

brier_tree <- read.csv("pred_prog_results_brier_cens.csv", header = FALSE)
brier_bart <- read.csv("../competingmethods/BART/pred_prog/barts_miss.csv") %>%
  mutate(method = "bart") %>%
  select(-brier) %>%
  pivot_longer(
    cols = c("brier_cens", "brier_cens_miss"),
    names_to = "type",
    values_to = "brier"
  ) %>%
  mutate(
    type = case_when(
      type == "brier_cens" ~ "all",
      type == "brier_cens_miss" ~ "subset"
    )
  )

brier_forest <- read.csv("../competingmethods/random_forest_pred_prog/forests_miss.csv") %>%
  mutate(method = "forest") %>%
  select(-brier) %>%
  pivot_longer(
    cols = c("brier_cens", "brier_cens_miss"),
    names_to = "type",
    values_to = "brier"
  ) %>%
  mutate(
    type = case_when(
      type == "brier_cens" ~ "all",
      type == "brier_cens_miss" ~ "subset"
    )
  )

brier_cox <- read.csv("pred_prog_results_brier_cens_miss_cox_model.csv", header = FALSE) %>%
  mutate(method = "cox") %>%
  pivot_longer(
    -method,
    names_to = "t",
    values_to = "brier"
  ) %>%
  mutate(
    t = case_when(
      t == "V1" ~ .55,
      t == "V2" ~ .81,
      t == "V3" ~ 1.13,
      t == "V4" ~ 1.84,
      t == "V5" ~ 2.50
    ),
    type = "subset"
  )

brier <- brier_tree %>%
  pivot_longer(
    everything(),
    names_to = "t",
    values_to = "brier_cens_miss"
  ) %>%
  mutate(
    t = case_when(
      t == "V1" ~ .55,
      t == "V2" ~ .81,
      t == "V3" ~ 1.13,
      t == "V4" ~ 1.84,
      t == "V5" ~ 2.50
    ),
    method = "tree",
    type = "subset",
    brier = brier_cens_miss
  ) %>%
  bind_rows(brier_bart, brier_forest, brier_cox) %>%
  mutate(
    percentile = case_when(
      t == 0.55 ~ "10%",
      t == 0.81 ~ "25%",
      t == 1.13 ~ "50%",
      t == 1.84 ~ "75%",
      t == 2.50 ~ "90%"
    ),
    Method = case_when(
      method == "bart" ~ "BART",
      method == "cox" ~ "Cox",
      method == "forest" ~ "Forest",
      method == "tree" ~ "THM"
    ),
    Method = factor(Method, c("Cox", "BART", "Forest", "THM"))
  )


p <- ggplot(
  filter(brier, type == "subset"),
  aes(percentile, brier, fill = Method)
) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_grey(start = 0.3, end = 0.9) +
  labs(y = "Brier Score", x = "Survival Percentile")

ggsave("pred_prog_brier.jpg", plot = p, width = 7, height = 3)

# calculate average censoring
files <- dir_ls("../competingmethods/data/")
files <- files[grepl("pred_prog_2[[:digit:]]{2}.csv", files)]
get_cens_rate <- function(file) {
  data <- read.csv(file)
  mean(data$Y_2 == 0)
}
cens_rate <- map_dbl(
  files,
  get_cens_rate
)
mean(cens_rate)

