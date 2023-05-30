library(dplyr)
library(tidyr)
library(ggplot2)

bias_tree <- read.csv("pred_prog_results_bias.csv", header = FALSE)
coverage_tree <- read.csv("pred_prog_results_coverage.csv", header = FALSE)
brier_tree <- read.csv("pred_prog_results_brier_cens.csv", header = FALSE)

bart <- read.csv("../competingmethods/BART/pred_prog/barts.csv") %>%
  mutate(method = "bart")
brier_bart <- read.csv("../competingmethods/BART/pred_prog/barts_miss.csv") %>%
  mutate(method = "bart")

forest <- read.csv("../competingmethods/random_forest_pred_prog/forests.csv") %>%
  mutate(method = "forest")
brier_forest <- read.csv("../competingmethods/random_forest_pred_prog/forests_miss.csv") %>%
  mutate(method = "forest")

brier_cox <- read.csv("pred_prog_results_brier_cens_cox_model.csv", header = FALSE) %>%
  mutate(method = "cox") %>%
  pivot_longer(
    -method,
    names_to = "t",
    values_to = "brier"
  ) %>%
  mutate(
    t = case_when(
      t == "V1" ~ .53,
      t == "V2" ~ .80,
      t == "V3" ~ 1.11,
      t == "V4" ~ 1.86,
      t == "V5" ~ 2.53
    )
  )

bias_cox <- read.csv("pred_prog_results_bias_cox_model.csv", header = FALSE) %>%
  mutate(method = "cox", obs = 1:n()) %>%
  pivot_longer(
    -all_of(c("method", "obs")),
    names_to = "t",
    values_to = "bias"
  ) %>%
  mutate(
    t = case_when(
      t == "V1" ~ .53,
      t == "V2" ~ .80,
      t == "V3" ~ 1.11,
      t == "V4" ~ 1.86,
      t == "V5" ~ 2.53
    )
  ) %>%
  group_by(t, method, obs) %>%
  summarize(bias = mean(bias), .groups = "drop") %>%
  select(-obs)


# bias
bias_bart <- select(bart, t, obs, bias, method) %>%
  group_by(t, obs, method) %>%
  summarize(bias = mean(bias), .groups = "drop") %>%
  select(-obs)
bias_forest <- select(forest, t, obs, bias, method) %>%
  group_by(t, obs, method) %>%
  summarize(bias = mean(bias), .groups = "drop") %>%
  select(-obs)
bias <- bias_tree %>%
  pivot_longer(
    everything(),
    names_to = "t",
    values_to = "bias"
  ) %>%
  mutate(
    t = case_when(
      t == "V1" ~ .53,
      t == "V2" ~ .80,
      t == "V3" ~ 1.11,
      t == "V4" ~ 1.86,
      t == "V5" ~ 2.53
    ),
    method = "tree"
  ) %>%
  bind_rows(bias_bart, bias_forest, bias_cox)
ggplot(bias, aes(factor(t), bias, fill = method)) +
  geom_boxplot()

brier <- brier_tree %>%
  pivot_longer(
    everything(),
    names_to = "t",
    values_to = "brier_cens_miss"
  ) %>%
  mutate(
    t = case_when(
      t == "V1" ~ .53,
      t == "V2" ~ .80,
      t == "V3" ~ 1.11,
      t == "V4" ~ 1.86,
      t == "V5" ~ 2.53
    ),
    method = "tree"
  ) %>%
  bind_rows(brier_bart, brier_forest, brier_cox)


ggplot(brier, aes(factor(t), brier_cens_miss, fill = method)) +
  geom_boxplot()


# coverage
coverage_bart <- select(bart, t, obs, coverage, method) %>%
  group_by(t, obs, method) %>%
  summarize(coverage = mean(coverage), .groups = "drop") %>%
  select(-obs)
coverage <- coverage_tree %>%
  pivot_longer(
    everything(),
    names_to = "t",
    values_to = "coverage"
  ) %>%
  mutate(
    t = case_when(
      t == "V1" ~ .55,
      t == "V2" ~ .81,
      t == "V3" ~ 1.13,
      t == "V4" ~ 1.84,
      t == "V5" ~ 2.50
    ),
    method = "tree"
  ) %>%
  bind_rows(coverage_bart)
ggplot(coverage, aes(factor(t), coverage, fill = method)) +
  geom_boxplot()

