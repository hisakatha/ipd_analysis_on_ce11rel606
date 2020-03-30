#!/usr/bin/env bash
OUTDIR=ame.low_ipd.ab.out
CONTROL=../../deep_kinetics_region.ab.bed.merged.sorted.slop20.c_elegans.bed.fa
INPUT=../../extreme_ipd.low.ab.coverage25.slop20.fullLength.c_elegans.gff.fa
MOTIF=motifs_low_ipd.meme
~/memesuite_5_1_1_install/bin/ame -oc $OUTDIR --control $CONTROL $INPUT $MOTIF
