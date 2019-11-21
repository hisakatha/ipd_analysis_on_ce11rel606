#!/usr/bin/env bash
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -N memechip_l_normBy_cd
#$ -pe smp 8
#$ -l hostname=ax*

echo "Adenine modification l_normBy_cd" > description
~/memesuite/bin/meme-chip -oc . -ccut 100 -fdesc description -order 1 -db ~/memesuite/db/motif_databases/WORM/uniprobe_worm.meme -meme-mod anr -meme-minw 4 -meme-maxw 30 -meme-nmotifs 8 -meme-searchsize 100000 -meme-p 8 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 ../../modifications.l_normBy_cd.score20_ipdratio2_coverage25.A.slop20.fullLength.gff.fa
