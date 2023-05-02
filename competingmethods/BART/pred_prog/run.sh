#!/bin/bash

module load R
nohup Rscript run.R > run.out 2> run.err &
