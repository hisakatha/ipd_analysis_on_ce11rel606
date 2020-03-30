#!/usr/bin/env bash
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -N memechip
#$ -pe smp 10
#$ -l hostname=ax*

MEMEPREFIX="/home/hisakatha/memesuite_5_1_1_install"
SEED=1234
echo "High IPD in C. elegans, sample ab (VC2010+OP50(WGA))" > description
$MEMEPREFIX/bin/meme-chip -oc result -seed $SEED -ccut 100 -fdesc description -order 1 -db $MEMEPREFIX/motif_databases/WORM/uniprobe_worm.meme -meme-mod anr -meme-minw 4 -meme-maxw 30 -meme-nmotifs 8 -meme-searchsize 100000 -meme-p 10 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 -neg ../../deep_kinetics_region.ab.bed.merged.sorted.slop20.c_elegans.bed.fa ../../extreme_ipd.high.ab.coverage25.slop20.fullLength.c_elegans.gff.fa
