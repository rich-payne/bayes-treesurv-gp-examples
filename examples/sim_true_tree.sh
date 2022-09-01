#!/bin/bash

#NOTE: SGE recognizes lines that start with #$ even though bash sees the # as a comment.
#      -N tells SGE the name of the job.
#      -o tells SGE the name of the output file.
#      -e tells SGE the name of the error file.
#      -cwd tells SGE to execute in the current working directory (cwd).

#SGE options
#$ -N sim_true_tree
#$ -o sim_true_tree.out
#$ -e sim_true_tree.err
#$ -cwd

module load matlab
sh sim_true_tree_gen_grid.sh
sh sim_true_tree_grid.sh
echo "finished"
