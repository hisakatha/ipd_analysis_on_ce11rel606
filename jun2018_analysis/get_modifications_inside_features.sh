#!/bin/bash

MOD=$1
MOD_BASE=$(basename $MOD .gff)

# Merged set of intron_outside_exon
INTRON_OUTSIDE_EXON="/glusterfs/hisakatha/wormbase/c_elegans.PRJNA13758.WS267.annotations.gff3.gz.ucsc_chr.intron_outside_exon"

# Non-merged introns
INTRON="/glusterfs/hisakatha/wormbase/c_elegans.PRJNA13758.WS267.annotations.gff3.gz.ucsc_chr.intron.sorted"

#INTRON_OUTSIDE_EXON_BASE=$(basename $INTRON_OUTSIDE_EXON .bed)
MOD_INSIDE_INTRON=$MOD_BASE.inside.intron_outside_exon

# Setting of strand for bedtools behavior
# "-s" turns on requiring the same strand overlaps
STRAND=""
#STRAND="-s"

NAME2GENE="/glusterfs/hisakatha/wormbase/c_elegans.PRJNA13758.WS267.xrefs.name2gene.sh"
#NAME2SYMBOL="/glusterfs/hisakatha/wormbase/c_elegans.PRJNA13758.WS267.xrefs.name2symbol.sh"

bedtools intersect $STRAND -a $MOD -b $INTRON_OUTSIDE_EXON > $MOD_INSIDE_INTRON
bedtools intersect $STRAND -a $MOD_INSIDE_INTRON -b $INTRON -wb | cut -f 18 | sed -E "s/^Parent=Transcript:([-_.0-9a-zA-Z]+).*/\1/" | sort -V | uniq > $MOD_INSIDE_INTRON.names
$NAME2GENE $MOD_INSIDE_INTRON.names | sort -V | uniq > $MOD_INSIDE_INTRON.names.gene
#$NAME2SYMBOL $MOD_INSIDE_INTRON.names | sort -V | uniq > $MOD_INSIDE_INTRON.names.symbol

TTS_DOWN_SORTED="c_elegans.PRJNA13758.WS267.annotations.tts_down500.sorted.bed"
TTS_DOWN_MERGED="c_elegans.PRJNA13758.WS267.annotations.tts_down500.sorted.merged.bed"
#TBASE=$(basename $TTS_M .bed)
MOD_INSIDE_TTS_DOWN=$MOD_BASE.inside.tts_down500

bedtools intersect $STRAND -a $MOD -b $TTS_DOWN_MERGED > $MOD_INSIDE_TTS_DOWN
bedtools intersect $STRAND -a $MOD_INSIDE_TTS_DOWN -b $TTS_DOWN_SORTED -wb | cut -f 13 | sed -E "s/^Parent=Transcript:([-_.0-9a-zA-Z]+).*/\1/" | sort -V | uniq > $MOD_INSIDE_TTS_DOWN.names
$NAME2GENE $MOD_INSIDE_TTS_DOWN.names | sort -V | uniq > $MOD_INSIDE_TTS_DOWN.names.gene
#$NAME2SYMBOL $MOD_INSIDE_TTS_DOWN.names | sort -V | uniq > $MOD_INSIDE_TTS_DOWN.names.symbol

