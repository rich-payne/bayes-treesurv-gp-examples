library(targets)
library(tarchetypes)

tar_option_set(
  packages = c("randomForestSRC", "dplyr", "survival", "tidyr", "purrr"),
  retrieval = "worker",
  storage = "worker"
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "sge")
options(clustermq.template = "clustermq.tmpl")

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

filenames <- sprintf("../data/pred_prog_%d.csv", 201:250)
filenames_indep <- sprintf("../data/pred_prog_%d.csv", 101:150)
data_files <- tar_map(
  list(filenames = filenames, seed = 201:250),
  names = "seed",
  tar_target(file, filenames, format = "file")
)
data_files_indep <- tar_map(
  list(filenames = filenames_indep, seed = 101:150),
  names = "seed",
  tar_target(file, filenames, format = "file")
)
forest_targets <- tar_map(
  values = list(
    file = rlang::syms(paste0("file_", 201:250)),
    file_indep = rlang::syms(paste0("file_", 101:150)),
    seed = 201:250,
    sim = 1:50
  ),
  names = "seed",
  tar_target(forest, run_forest(file, file_indep, data_file_ref, sim))
  # tar_target(brier_cens_miss, get_brier_miss(forest, missing_data, sim))
)

write_file <- function(x, file) {
  write.csv(x, file, row.names = FALSE)
  file
}

list(
  data_files,
  data_files_indep,
  tar_target(data_file_ref, "../data/pred_prog_12345.csv", format = "file"),
  tar_target(missing_data_file, "../../examples/pred_prog_results_missing_detail.csv", format = "file"),
  tar_target(missing_data, read.csv(missing_data_file)),
  forest_targets,
  tar_combine(forests, forest_targets$forest, command = bind_rows(!!!.x)),
  tar_target(forests_miss, add_brier_miss(forests, missing_data)),
  tar_target(forests_file, write_file(forests, "forests.csv"), format = "file"),
  tar_target(forests_file_miss, write_file(forests_miss, "forests_miss.csv"), format = "file")
)
