#!/bin/bash

experiment_name="permval"
#for marginal in gauss ng 
#do
#  for cordiff in gap order rand
#  do
    #for i in 1 2 3 4
    #do
      #for s in $(seq 1  200)
      for s in 107 112 114 115 117 118 127 131 134 137 142 148 157 168 182 195 28 63 71 77 85 92
      do
        echo "Permutation test of instance "$i", perm batch = "$s
        command="$Rscript ./exp/perm.sh $s"
        sbatch --time=01:00:00 -p RM-shared -N 1 -n 1 -J $experiment_name $command
      done
    #done
# done
#done
