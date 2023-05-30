#!/bin/bash

module load R # Uncomment if R is an environment module.
nohup Rscript run.R > run.out 2> run.err &
