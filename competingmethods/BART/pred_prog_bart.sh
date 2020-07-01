#!/bin/bash

#NOTE: SGE recognizes lines that start with #$ even though bash sees the # as a comment.
#      -N tells SGE the name of the job.
#      -o tells SGE the name of the output file.
#      -e tells SGE the name of the error file.
#      -cwd tells SGE to execute in the current working directory (cwd).

#SGE options
#$ -N validation
#$ -o validation.out
#$ -e validation.err
#$ -cwd
# #$ -pe smp 4
# #$ -R y

module load R-qualified
nohup Rscript --no-restore pred_prog_bart.R > pred_prog_bart.out 2> pred_prog_bart.err &
