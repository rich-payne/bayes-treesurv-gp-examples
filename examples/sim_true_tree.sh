#!/bin/bash

module load qarray
sh sim_true_tree_gen_grid.sh
nohup qarray -f sim_true_tree_grid.sh > sim_true_tree.out 2> sim_true_tree.err &
