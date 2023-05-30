library(targets)
library(tarchetypes)

# Set target options:
tar_option_set(
  packages = c("BART", "dplyr", "randomForestSRC", "tibble", "tidyr"),
  retrieval = "worker",
  storage = "worker",
  error = "abridge"
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "sge")
options(clustermq.template = "clustermq.tmpl")

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

kfold_targs <- tar_map(
  values = list(fold = 1:10),
  tar_target(forest_fit, fit_forest(data, k_fold = fold)),
  tar_target(bart_fit, fit_bart(data, k_fold = fold)),
  tar_target(
    brier_forest,
    get_brier_forest(
      forest_fit,
      times = c(617, 1181, 1788, 2691, 3608),
      data = data,
      fold = fold,
      missing_data = missing_data
    )
  ),
  tar_target(
    brier_bart,
    get_brier_bart(
      bart_fit,
      times = c(617, 1181, 1788, 2691, 3608),
      data = data,
      fold = fold,
      missing_data = missing_data
    )
  )
)

list(
  tar_target(pbc_file, "../../examples/pbc_kfold.csv", format = "file"),
  tar_target(data, read.csv(pbc_file, stringsAsFactors = TRUE)),
  tar_target(missing_data_file, "../../examples/pbc_analysis_results_missing_detail.csv", format = "file"),
  tar_target(missing_data, read.csv(missing_data_file)),
  tar_target(brier_tree_file, "../../examples/pbc_analysis_results_brier_cens.csv", format = "file"),
  tar_target(brier_tree, read.csv(brier_tree_file, header = FALSE)),
  tar_target(forest_fit, fit_forest(data)),
  tar_target(bart_fit, fit_bart(data)),
  kfold_targs,
  tar_combine(brier_forest, kfold_targs$brier_forest, command = bind_rows(!!!.x)),
  tar_combine(brier_bart, kfold_targs$brier_bart, command = bind_rows(!!!.x)),
  tar_render(report, "pbc.Rmd")
)
