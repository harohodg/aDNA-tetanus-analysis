#!/usr/bin/env bash


sucessful=$(grep "DONE" *_fastp*.out | wc -l)
running=$(sq -h --format="%j" --state='RUNNING' | grep "fastp" | wc -l)
pending=$(sq -h --format="%j" --state='PENDING' | grep "fastp" | wc -l)
out_of_mem=$(grep "killed by the cgroup out-of-memory" *_fastp*.out | wc -l)
out_of_time=$(grep "DUE TO TIME LIMIT" *_fastp*.out | wc -l)
canceled=$(grep "CANCELLED AT.*" *_fastp*.out | grep -v "DUE TO TIME LIMIT" | wc -l)
echo "Successful : $sucessful "
echo "Running    : $running"
echo "Pending    : $pending"
echo "Out of memory : $out_of_mem "
echo "Out of Time : $out_of_time"
echo "Canceled : $canceled"

echo "Total $(expr $sucessful + $running + $pending + $out_of_mem + $out_of_time + $canceled)"
