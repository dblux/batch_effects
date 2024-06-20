#!/bin/bash

# Command-line arguments
# $1: size 
# $2: outdir


logfile="$2log-n$1.txt"
timefile="$2time-n$1.txt"
memfile="$2memory-n$1.txt"
echo "LOG FILE" > $logfile
echo "metric n user system elapsed user.child sys.child" > $timefile 
echo "metric n maxsize" > $memfile 
echo "Created file: $logfile"
echo "Created file: $timefile"
echo "Created file: $memfile"

# sizes=($(seq 2000 2000 10000))
# metrics=(HVP gPCA PVCA CMS kBET LISI)
metrics=(HVP)

for metric in ${metrics[@]}
do
  echo "Benchmarking: $metric (n = $1)"
  /usr/bin/time -f "%M" \
    Rscript notebooks/benchmark/benchmark.R $metric $1 $2 \
    2>&1 >>$timefile | \
    tee -a $logfile | \
    awk -v m="$metric" -v s="$1" 'END { print m,s,$1 }' \
    >>$memfile
done
