#!/bin/bash
rm pbc_analysis_grid.sh
for i in {1..10};
do
  echo "./pbc_analysisstandaloneApplication/run_pbc_analysis.sh <insert-path-to-matlab> 'pbc_kfold.csv' $i $i 2" >> pbc_analysis_grid.sh
done
