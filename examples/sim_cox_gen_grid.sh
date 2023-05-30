#!/bin/bash
rm sim_cox_grid.sh
for i in {1..50};
do
  echo "./sim_coxstandaloneApplication/run_sim_cox.sh /lrlhps/apps/bioinfo/matlab/matlab_2022a $((100 + $i)) $i 10000 10" >> sim_cox_grid.sh
done
