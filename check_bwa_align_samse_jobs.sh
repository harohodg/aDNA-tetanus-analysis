#!/usr/bin/env bash


sucessful=$(grep "Done" *_bwa_align_samse*.out | wc -l)
running=$(sq -h --format="%j" --state='RUNNING' | grep "bwa_align_samse" | wc -l)
pending=$(sq -h --format="%j" --state='PENDING' | grep "bwa_align_samse" | wc -l)
out_of_mem=$(grep "killed by the cgroup out-of-memory" *_bwa_align_samse*.out | wc -l)
out_of_time=$(grep "DUE TO TIME LIMIT" *_bwa_align_samse*.out | wc -l)
canceled=$(grep "CANCELLED AT.*" *_bwa_align_samse*.out | grep -v "DUE TO TIME LIMIT" | wc -l)
echo "Successful : $sucessful "
echo "Running    : $running"
echo "Pending    : $pending"
echo "Out of memory : $out_of_mem "
echo "Out of Time : $out_of_time"
echo "Canceled : $canceled"

echo "Total $(expr $sucessful + $running + $pending + $out_of_mem + $out_of_time + $canceled)"
