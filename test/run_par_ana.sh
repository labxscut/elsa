#!/bin/bash
#par_ana to 30 cores; 5M comparison per core
pf=cell_umi_cts
echo "par_ana $pf.txt $pf.lsa 'lsa_compute %s %s -e $pf.txt -n rnz -b 0 -s 2789 -r 1 -p theo' $PWD >$pf.par.log 2>&1"
cmd="par_ana $pf.txt $pf.lsa 'lsa_compute %s %s -e $pf.txt -n rnz -b 0 -s 2789 -r 1 -p theo' $PWD >$pf.par.log 2>&1"
eval "$cmd"
