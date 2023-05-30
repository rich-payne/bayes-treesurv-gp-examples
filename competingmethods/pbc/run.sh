#!/bin/bash

module load R-qualified
nohup Rscript run.R > run.out 2> run.err &
