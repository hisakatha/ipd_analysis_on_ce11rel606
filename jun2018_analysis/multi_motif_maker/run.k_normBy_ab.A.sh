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
#MINSCORE=0

GFF="../modifications.k_normBy_ab.score20_ipdratio2_coverage25.A.withContext.gff"
GFF2="$(basename $GFF .gff).for_motifmaker.gff"
# append coverage=min(native_coverage, control_coverage) for MotifMaker
(cat $GFF | grep "^#"; cat $GFF | grep -v "^#" | awk -v OFS="\t" \
'{split($9, arr, ";"); if(substr(arr[2], 1, 15) != "native_coverage"){exit 11}; if(substr(arr[3], 1, 16) != "control_coverage"){exit 12};'\
' split(arr[2], nat_cov_arr, "="); split(arr[3], ctl_cov_arr, "=");'\
' if(nat_cov_arr[2] < ctl_cov_arr[2]){min=nat_cov_arr[2]}else{min=ctl_cov_arr[2]}; $9 = $9 ";coverage=" min; print $0}') > $GFF2

java -Xmx128G -jar $JAR find -t $THREAD -f $REF -g $GFF2 -o "$(basename $GFF2).motifs.csv"

