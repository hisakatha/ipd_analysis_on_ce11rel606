#!/usr/bin/env bash
# gstart, gend, gstrand: start position, end position, and strand of an input genome sequence
# mstart, mend, mstrand: start position, end position, and strand of a detected motif

# output: chr, motif start in bed format (0-based), strand
# For both + and - strand, motif end = motif start + motif length
cat $1 | sed -E 's/^([-\.\+=|_a-zA-Z0-9]+)(:([0-9]+)-([0-9]+)(\([-\+\.]\))?)?\t([0-9]+)\t([0-9]+).+([-\+])$/\1\t\3\t\4\t\5\t\6\t\7\t\8/' |
 awk -v FS="\t" '{chr=$1; gstart=$2; gend=$3; gstrand=$4; mstart=$5; mend=$6; mstrand=$7; if(gstrand != "(-)"){print chr,gstart+mstart,mstrand}else{print chr,gend-mend,((mstrand=="+")?"-":"+")}}' |
 sort | uniq
