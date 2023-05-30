#!/bin/bash

module load qarray
sh sim_cox_gen_grid.sh
nohup qarray -f sim_cox_grid.sh > sim_cox.out 2> sim_cox.err &
