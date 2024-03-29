These files contain the code to reproduce the examples in the manuscript "A Bayesian Survival Treed Hazards Model Using Latent Gaussian Processes" accepted in *Biometrics*.

The bayes-treesurv-gp directory contains the code to run the treed hazards model.

Code in the examples/ directory contains the code to apply the treed hazards model (THM) and the Cox model on the primary biliary cirrhosis dataset (pbc) and the predictive/prognostic simulations.  The files include:

* Tutorial.m: A tutorial on how to run the THM code.
* gen_data_pred_prog.m: Generate data for the predictive/prognostic simulations.
* get_brier_score.m, get_brier_score_cens.m: Functions for calculating the brier score with/without censored observations.
* get_x_pbc.m: Create the design matrix for the pbc analysis.
* pbc.R: Clean and write the pbc datasets to file.
* pbc.m: Run the THM analysis on the pbc dataset.
* pbc_analysis.m: Run the THM on a single fold of the pbc data.
* pbc_analysis.sh: Run the pbc analysis on each fold.
* pbc_analysis_compiler.m: Code to compile the code to run in parallel in a high performance computing environment.
* pbc_analysis_cox.m: Code to run the pbc data with a Cox model.
* pbc_analysis_gen_grid.sh: Generate the pbd analysis grid for k-fold validation.
* pbc_analysis_results.m: Summarize results from the pbc analyses.
* pred_prog_cox.m: Run the Cox analysis on the predictive/prognostic simulations.
* pred_prog_results.R: Create the Brier score plots for the predictive/prognostic simulations.
* pred_prog_results_single_tree.m: Obtain the figure from the first predictive/prognostic simulation.
* run_mcmc.m: A helper function to run the MCMC for the THM.
* sim_pred_prog.m: A function to run a single predictive/prognostic simulation.
* sim_pred_prog.sh: Run the predictive/prognostic simulations on the cluster.
* sim_pred_prog_compiler.m: Code to compile the code to run in parallel in a high performance computing environment.
* sim_pred_prog_gen_grid.sh: Generate the simulation grid.
* sim_pred_prog_results.m: Summarize the results of the predictive/prognostic simulations.
* write_pred_prog_data.m: Write the datasets for the predictive/prognostic simulations.

In the competingmethods/ directory, the pbc analyses and predictive/prognostic simulations are analyzed with BART and random forest methods.  The pbc/ directory contains a targets pipeline which analyzes the pbc k-folds for both BART and random forest.  The BART/ and random_forest_pred_prog/ directories run BART and random forest on the predictive/prognostic simulations, respectively, also utilizing targets pipelines.
