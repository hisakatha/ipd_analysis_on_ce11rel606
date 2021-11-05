#!/usr/bin/env bash
out="count_occ_simple.csv"
header_fragment="motif,sample_id"
cat GATC/count_occ_simple_per_motif.csv | head -n1 > $out
cat */count_occ_simple_per_motif.csv | grep -v "$header_fragment" >> $out
