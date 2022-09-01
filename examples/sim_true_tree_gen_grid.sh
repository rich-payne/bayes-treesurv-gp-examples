#!/bin/bash
rm sim_true_tree_grid.sh
for i in {1..50};
do
  echo "matlab -nodisplay -nodesktop -r \"sim_true_tree($i, $i); quit;\"" >> sim_true_tree_grid.sh
done
