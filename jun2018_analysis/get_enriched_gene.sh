#!/bin/bash

DATABASE_INPUTS=$(ls /glusterfs/hisakatha/wormbase/c_elegans.PRJNA13758.WS267.annotations.gff3.gz.ucsc_chr.{exon,intron,ncrna,tss_up500,tss_down500,tts_up500,tts_down500}.sorted)
DATABASE=get_enriched_gene.feature_database
#echo "chr,feature,start,end,strand,transcript,gene" > $DATABASE
echo -e "chr\tfeature\tstart\tend\tstrand\ttranscript\tgene" > $DATABASE
NAME2GENE="/glusterfs/hisakatha/wormbase/c_elegans.PRJNA13758.WS267.xrefs.name2gene.sh"
#paste -d, <(cat $DATABASE_INPUTS | cut -f 1,3,4,5,7 --output-delimiter=, ) \
#    <(cat $DATABASE_INPUTS | sed -E "s/.*Transcript:([-_\.0-9a-zA-Z]+).*/\1/") \
#    <(cat $DATABASE_INPUTS | sed -E "s/.*Transcript:([-_\.0-9a-zA-Z]+).*/\1/" | $NAME2GENE ) >> $DATABASE

paste <(cat $DATABASE_INPUTS | cut -f 1,3,4,5,7 ) \
    <(cat $DATABASE_INPUTS | sed -E "s/.*Transcript:([-_\.0-9a-zA-Z]+).*/\1/") \
    <(cat $DATABASE_INPUTS | sed -E "s/.*Transcript:([-_\.0-9a-zA-Z]+).*/\1/" | $NAME2GENE ) >> $DATABASE
