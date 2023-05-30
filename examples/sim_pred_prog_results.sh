#!/bin/bash

module load matlab
nohup matlab -nodisplay -r "sim_pred_prog_results; exit" > sim_pred_prog_results.out 2> sim_pred_prog_results.err &
