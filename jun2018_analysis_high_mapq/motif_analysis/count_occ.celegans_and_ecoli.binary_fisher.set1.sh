#!/usr/bin/env bash

motifs_high_ipd="ACGCRTG ATCAGCTG GGN_4 RGTA AAAABT AGAGAGTA AGCTATAT AGGCAGGC ATGGGAYA ATTGTTAC CAAACTAC CAGYTG CCAATCAG CRACGAS DCGAGACC DGCTTC GAAGGATC GATATRGY GCGACCTA GCGCGCGC GGHGGY GTAGATCA GTATCGTA TGACGTCA TGGTGSA YGGAR"
n_motifs_high_ipd=26

motifs_low_ipd="ACATMTGG CTGDAR TACTGTAG AATMAATA AGACGCAG CGGYTTGA CTKCAA GCGCGTCA GTGTGTGY TACCCCKA TACCTTGA TGACGTCA"
n_motifs_low_ipd=12

motifs_select1="GAGG AGAA"

csv_prefix="count_occ.celegans_and_ecoli.binary_fisher.set1"
pre_csv="$csv_prefix.pre_csv"

motifs=$(echo $motifs_high_ipd $motifs_low_ipd $motifs_select1 | xargs -n1 | sort | uniq)
for motif in $motifs
do
    cd $motif
    for sample in ab cd k l PD2182 PD2182sequel
    do
        ../count_occ.celegans_and_ecoli.binary_fisher.sh $sample $n_motifs_high_ipd $n_motifs_low_ipd || exit 1
    done
    cd ../
done > $pre_csv

(cat $pre_csv | grep sample_id | head -n1; cat $pre_csv | grep -v "#" | grep -v sample_id) > $csv_prefix.csv
