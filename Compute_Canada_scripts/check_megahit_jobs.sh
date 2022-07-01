#!/usr/bin/env bash


sucessful=$(grep "ALL DONE" *_megahit*.out | wc -l)
running=$(sq -h --format="%j" --state='RUNNING' | grep "megahit" | wc -l)
pending=$(sq -h --format="%j" --state='PENDING' | grep "megahit" | wc -l)
#error=$(grep "Error occurs" *_megahit*.out | wc -l)
out_of_mem=$(grep "killed by the cgroup out-of-memory" *_megahit*.out | wc -l)
out_of_time=$(grep "DUE TO TIME LIMIT" *_megahit*.out | wc -l)
canceled=$(grep "CANCELLED AT.*" *_megahit*.out | grep -v "DUE TO TIME LIMIT" | wc -l)
echo "Successful : $sucessful "
echo "Running    : $running"
echo "Pending    : $pending"
#echo "Error      : $error"
echo "Out of memory : $out_of_mem "
echo "Out of Time : $out_of_time"
echo "Canceled : $canceled"

echo "Total $(expr $sucessful + $running + $pending + $out_of_mem + $out_of_time + $canceled)"
