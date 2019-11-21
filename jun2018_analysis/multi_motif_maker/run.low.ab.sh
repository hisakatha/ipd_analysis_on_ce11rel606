#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd
#$ -N motifmaker
#$ -j y
#$ -V
#$ -pe smp 24
#$ -l hostname=ax*
#$ -l mem_free=128G

REF="/glusterfs/hisakatha/ce11rel606/ce11rel606.fa"
THREAD=24
JAR="/home/hisakatha/MultiMotifMaker/artifacts/MultiMotifMaker.jar"
MINSCORE=0

#GFF="../modifications.l_normBy_cd.score20_ipdratio2_coverage25.A.withContext.gff"

GFF="../extreme_ipd.low.ab.coverage25.withContext.gff"
java -Xmx128G -jar $JAR find -t $THREAD -f $REF -g $GFF -m $MINSCORE -o "$(basename $GFF).motifs.csv"

