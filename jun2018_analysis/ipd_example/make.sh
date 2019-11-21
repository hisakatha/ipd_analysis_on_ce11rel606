#!/usr/bin/env bash
Rscript get_ipd_in_occ.R high_ipdratio_adenine2.occ
dna_seq=$(echo -e "chrIII\t12802497\t12802518\t.\t.\t-" | bedtools getfasta -s -fi /glusterfs/hisakatha/ce11rel606/ce11rel606.fa -bed - | grep -v ">")
Rscript plot_kinetics.stderr.R high_ipdratio_adenine2.occ.ipd $dna_seq
