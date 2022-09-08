#!/bin/bash
rm sim_true_tree_grid.sh
for i in {1..50};
do
  echo "./sim_true_treestandaloneApplication/run_sim_true_tree.sh /lrlhps/apps/bioinfo/matlab/matlab_2022a $i $i 10000 10" >> sim_true_tree_grid.sh
done
