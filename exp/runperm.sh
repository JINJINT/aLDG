#!/bin/bash

experiment_name="permval"
#for marginal in gauss ng 
#do
#  for cordiff in gap order rand
#  do
    #for i in 1 2 3 4
    #do
      for s in $(seq 1 2)
      do
        echo "Permutation test of instance "$i", perm batch = "$s
        command="$Rscript ./exp/perm.sh $s"
        sbatch --time=00:30:00 -p RM-shared -N 1 -n 1 -J $experiment_name $command
      done
    #done
# done
#done
