#!/usr/bin/env bash
#$ -S /bin/bash
#$ -V
#$ -j y
#$ -o ipd_summary_high_mapq.log
#$ -cwd
#$ -N ipdSummary
#$ -pe smp 24

export PATH="/bio/package/pacbio/smrtlink/smrtcmds/bin:$PATH"

fname="ipd_summary_high_mapq"
ipdSummary output/tasks/pbalign.tasks.pbalign-0/mapped.alignmentset.bam --reference ../ce11rel606.xml --gff $fname.gff --bigwig $fname.bw --csv_h5 $fname.h5 --numWorkers 24 --pvalue 0.001 --verbose --maxLength 3000000000000 --methylFraction --identify m6A,m4C --mapQvThreshold 128

# original command by pbsmrtpipe
#ipdSummary /glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_ab_pcr_vc2010_op50/output/tasks/pbcoretools.tasks.alignment_contig_scatter_basemods-1/chunk_alignmentset_0.alignmentset.xml --reference /glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/ce11rel606.xml --gff /glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_ab_pcr_vc2010_op50/output/tasks/kinetics_tools.tasks.ipd_summary-1/basemods.gff --bigwig /glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_ab_pcr_vc2010_op50/output/tasks/kinetics_tools.tasks.ipd_summary-1/basemods.bw --csv_h5 /glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_ab_pcr_vc2010_op50/output/tasks/kinetics_tools.tasks.ipd_summary-1/basemods.h5 --numWorkers 7 --pvalue 0.001 --alignmentSetRefWindows --verbose --maxLength 3000000000000 --methylFraction --identify m6A,m4C
