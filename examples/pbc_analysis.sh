#!/bin/bash

module load qarray
sh pbc_analysis_gen_grid.sh
nohup qarray -f pbc_analysis_grid.sh > pbc_analysis.out 2> pbc_analysis.err &
