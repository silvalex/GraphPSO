#!/bin/bash

datasetroot=~/workspace
datasets='wsc2008 wsc2009'
resultsroot=~/workspace/GraphPSOResults
for d in $datasets
do
    for t in $(find ${datasetroot}/$d/ -mindepth 1 -type d)
    do
      set=$(basename $t)
      result_dir=${set}Results
      mkdir -p $resultsroot/$result_dir
      for i in $(seq 0 29)
      do
	java -classpath ~/workspace/GraphPSO/bin pso.GraphPSO ${resultsroot}/${result_dir}/out${i}.stat ${t}/problem.xml ${t}/services-output.xml ${t}/taxonomy.xml ${i}
      done
    done
done





