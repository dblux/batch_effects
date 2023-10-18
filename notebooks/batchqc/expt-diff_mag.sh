#!/bin/bash

# CMDLINE ARGUMENTS
# $1: outdir for metric results


metrics=(RVP kBET LISI)
sizes=($(seq 0 1000 10000))


for size in ${sizes[@]}
do
  for metric in ${metrics[@]}
  do
    Rscript expt-diff_mag.R $metric $size $1
  done
done
