#!/bin/bash
rm sim_pred_prog_grid.sh
for i in {1..50};
do
  echo "./sim_pred_progstandaloneApplication/run_sim_pred_prog.sh <insert-path-to-matlab> $((200 + $i)) $i 10000 10" >> sim_pred_prog_grid.sh
done
