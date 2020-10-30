#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -j y
#$ -o pbsmrtpipe_log
#$ -cwd
#$ -N pbsmrtpipe

. /glusterfs/hisakatha/methylation/smrtpipe/setup.sh
mkdir -p output

pbsmrtpipe pipeline-id --entry "eid_subread:$(basename $(pwd)).subreads.xml" --entry eid_ref_dataset:../ce11rel606.xml --output-dir=./output --preset-xml=/bio/package/pacbio/smrtlink/userdata/config/preset.xml --preset-xml=../preset.ds_modification_detection.1.xml --disable-chunk-mode pbsmrtpipe.pipelines.ds_modification_detection

