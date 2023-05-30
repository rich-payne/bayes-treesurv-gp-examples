#!/bin/bash
rm pbc_analysis_grid.sh
for i in {1..10};
do
  echo "./pbc_analysisstandaloneApplication/run_pbc_analysis.sh /lrlhps/apps/bioinfo/matlab/matlab_2022a 'pbc_kfold.csv' $i $i 2" >> pbc_analysis_grid.sh
done
