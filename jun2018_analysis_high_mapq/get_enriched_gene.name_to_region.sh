#!/bin/bash
DATABASE="/glusterfs/hisakatha/wormbase/c_elegans.PRJNA13758.WS267.annotations.gff3.gz.ucsc_chr.gene"

# Make a list from a pipe or standard input
#GENES=$(cat)

GENES=$(cat get_enriched_gene.both.summarize.csv | cut -f1 -d, | uniq | tail -n +2)

for g in $GENES
do
    cat $DATABASE | grep -E "sequence_name=$g[;$]"
done

