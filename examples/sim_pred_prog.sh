#!/bin/bash

module load qarray
sh sim_pred_prog_gen_grid.sh
nohup qarray -f sim_pred_prog_grid.sh > sim_pred_prog.out 2> sim_pred_prog.err &
