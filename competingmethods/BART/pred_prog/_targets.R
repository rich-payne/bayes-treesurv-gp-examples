library(targets)
library(tarchetypes)
tar_option_set(
  packages = c("BART", "dplyr", "purrr", "tidyr", "survival"),
  # error = "abridge",
  storage = "worker",
  retrieval = "worker"
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "sge")
options(clustermq.template = "clustermq.tmpl")

tar_source()

filenames <- sprintf("../../data/pred_prog_%d.csv", 201:250)
filenames_indep <- sprintf("../../data/pred_prog_%d.csv", 101:150)
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
barts_targets <- tar_map(
  values = list(
    file = rlang::syms(paste0("file_", 201:250)),
    file_indep = rlang::syms(paste0("file_", 101:150)),
    seed = 201:250,
    sim = 1:50
  ),
  names = "seed",
  tar_target(bart, run_bart(file, file_indep, data_file_ref, sim))
)

write_file <- function(x, file) {
  write.csv(x, file, row.names = FALSE)
  file
}

# bind_rows_bart <- function(...) {
#   barts <- list(...)
#   stats <- map_dfr(barts, ~ .x$stats)
#   return(stats)
# }

list(
  data_files,
  data_files_indep,
  tar_target(data_file_ref, "../../data/pred_prog_12345.csv", format = "file"),
  tar_target(missing_data_file, "../../../examples/pred_prog_results_missing_detail.csv", format = "file"),
  tar_target(missing_data, read.csv(missing_data_file)),
  barts_targets,
  tar_combine(barts, barts_targets, command = bind_rows(!!!.x)),
  tar_target(barts_miss, add_brier_miss(barts, missing_data)),
  tar_target(barts_file, write_file(barts, "barts.csv"), format = "file"),
  tar_target(barts_file_miss, write_file(barts_miss, "barts_miss.csv"), format = "file")
  # tar_target(ocs, get_ocs(barts))
)
