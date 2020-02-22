#!/usr/bin/env bash
# Skip getting IPDs assuming make3.sh has been executed in advance
#Rscript get_ipd_in_occ.R high_ipdratio_adenine3.occ
dna_seq=$(echo -e "chrIII\t1319725\t1319746\t.\t.\t-" | bedtools getfasta -s -fi /glusterfs/hisakatha/ce11rel606/ce11rel606.fa -bed - | grep -v ">")
Rscript plot_kinetics.stderr.simple.R high_ipdratio_adenine3.occ.ipd $dna_seq
