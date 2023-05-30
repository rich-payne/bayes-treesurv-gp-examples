#!/bin/bash

module load matlab
nohup matlab -nodisplay -r "sim_true_tree_results; exit" > sim_true_tree_results.out 2> sim_true_tree_results.err &
