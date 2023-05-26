#!/bin/bash

# CMDLINE ARGUMENTS
# $1: outfile for profiling results
# $2: outdir for metric results


metrics=(RVP kBET LISI)
sizes=($(seq 2000 2000 10000))

echo "metric,samples,runtime,memory" > $1 

for size in ${sizes[@]}
do
  for metric in ${metrics[@]}
  do
    /usr/bin/time -f "%M" Rscript benchmark.R $metric $size $2 |& tee -a $1 
  done
done
