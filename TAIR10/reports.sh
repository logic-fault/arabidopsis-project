#!/bin/bash



seq_dir=raw_genome
report_dir=reports

for name in $(ls $seq_dir | grep \.fas)
do
   newcpgreport -sequence $seq_dir/$name -window 100 -shift 1 -minlen 200 -minoe 0.6 -minpc 50. -outfile $report_dir/$name.cpgreport
done
