#!/bin/bash

module load matlab
nohup matlab -nodisplay -r "sim_cox_results; exit" > sim_cox_results.out 2> sim_cox_results.err &
