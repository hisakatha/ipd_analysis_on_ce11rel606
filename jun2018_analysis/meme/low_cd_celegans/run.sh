#!/usr/bin/env bash
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -N memechip
#$ -pe smp 10
#$ -l hostname=ax*

MEMEPREFIX="/home/hisakatha/memesuite_5_1_1_install"
MEMEBIN="$MEMEPREFIX/bin"
MEMELIBEXEC="$MEMEPREFIX/libexec/meme-5.1.1"

SEED=1234
echo "Low IPD in C. elegans, sample cd (VC2010(WGA))" > description
CONTROL=../../deep_kinetics_region.cd.bed.merged.sorted.slop20.c_elegans.bed.fa

$MEMELIBEXEC/fasta-get-markov -m 1 -dna $CONTROL control_background
$MEMEBIN/meme-chip -oc result -seed $SEED -ccut 100 -fdesc description -bfile control_background -db $MEMEPREFIX/motif_databases/WORM/uniprobe_worm.meme -meme-mod anr -meme-minw 4 -meme-maxw 30 -meme-nmotifs 8 -meme-searchsize 100000 -meme-p 10 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 ../../extreme_ipd.low.cd.coverage25.slop20.fullLength.c_elegans.gff.fa
