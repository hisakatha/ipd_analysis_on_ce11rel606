#!/usr/bin/env bash
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -N meme
#$ -pe smp 22
#$ -l hostname=ax*

MEMEPREFIX="/home/hisakatha/memesuite_5_1_1_install"
MEMEBIN="$MEMEPREFIX/bin"
MEMELIBEXEC="$MEMEPREFIX/libexec/meme-5.1.1"

#SEED=1234
#echo "High IPD in C. elegans, sample ab (VC2010+OP50(WGA))" > description
#CONTROL=../../deep_kinetics_region.ab.bed.merged.sorted.slop20.c_elegans.bed.fa

$MEMEBIN/meme ../result/extreme_ipd.high.ab.coverage25.slop20.fullLength.c_elegans.gff.fa -oc meme_out -mod zoops -nmotifs 8 -minw 4 -maxw 30 -bfile ../result/control_background -dna -searchsize 1000000 -revcomp -nostatus -p 22
