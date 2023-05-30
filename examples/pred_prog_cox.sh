#!/bin/bash

module load matlab
nohup matlab -nodisplay -r "pred_prog_cox; exit" > pred_prog_cox.out 2> pred_prog_cox.err &
