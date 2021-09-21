#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -j y
#$ -o pbsmrtpipe_log
#$ -cwd
#$ -N pbsmrtpipe

export PATH="/bio/package/pacbio/smrtlink/smrtcmds/bin:$PATH"
mkdir -p output

pbsmrtpipe pipeline-id --entry "eid_subread:$(basename $(pwd)).subreads.xml" --entry eid_ref_dataset:../ce11rel606.xml --output-dir=./output --preset-xml=../preset.xml --preset-xml=../preset.ds_modification_detection.1.xml --disable-chunk-mode pbsmrtpipe.pipelines.ds_modification_detection

