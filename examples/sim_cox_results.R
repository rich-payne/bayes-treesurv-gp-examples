library(dplyr)
library(tidyr)
library(ggplot2)

bias_tree <- read.csv("cox_results_bias.csv", header = FALSE)
coverage_tree <- read.csv("cox_results_coverage.csv", header = FALSE)
brier_tree <- read.csv("cox_results_brier.csv", header = FALSE)

bart <- read.csv("../competingmethods/BART/cox/barts.csv") %>%
  mutate(method = "bart")
brier_bart <- distinct(bart, t, brier, method) %>%
  bind_rows(
    distinct(bart, t, brier_true, method) %>%
      mutate(method = "true") %>%
      rename(brier = brier_true)
  )

# bias
bias_bart <- select(bart, t, obs, bias, method) %>%
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
      t == "V1" ~ .11,
      t == "V2" ~ .36,
      t == "V3" ~ 1.24,
      t == "V4" ~ 4.08,
      t == "V5" ~ 10.36
    ),
    method = "tree"
  ) %>%
  bind_rows(bias_bart)
ggplot(bias, aes(factor(t), bias, fill = method)) +
  geom_boxplot()


brier <- brier_tree %>%
  pivot_longer(
    everything(),
    names_to = "t",
    values_to = "brier"
  ) %>%
  mutate(
    t = case_when(
      t == "V1" ~ .11,
      t == "V2" ~ .36,
      t == "V3" ~ 1.24,
      t == "V4" ~ 4.08,
      t == "V5" ~ 10.36
    ),
    method = "tree"
  ) %>%
  bind_rows(brier_bart)


ggplot(brier, aes(factor(t), brier, fill = method)) +
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
      t == "V1" ~ .11,
      t == "V2" ~ .36,
      t == "V3" ~ 1.24,
      t == "V4" ~ 4.08,
      t == "V5" ~ 10.36
    ),
    method = "tree"
  ) %>%
  bind_rows(coverage_bart)
ggplot(coverage, aes(factor(t), coverage, fill = method)) +
  geom_boxplot()

