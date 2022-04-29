#!/bin/bash/

i=1
while IFS="," read -r rec_column1 rec_column2 rec_column3 rec_column4
do
  echo "run_$i"
  check=$(expr $i % 10) 
  Rscript --vanilla Simulations.R $rec_column1 $rec_column2 $rec_column3 $rec_column4 &

  if [ "$check" -eq 0 ]
  then
	  sleep 23m
  fi
  ((i=i+1))
done < sim_params_missing.csv

wait
echo "All Done"

