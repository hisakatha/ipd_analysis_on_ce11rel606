#!/usr/bin/env bash
ln -s output/tasks/pbalign.tasks.pbalign-0/mapped.alignmentset.bam mapped.alignmentset.merged.bam
mkdir -p output/tasks/kinetics_tools.tasks.gather_kinetics_h5-1
ln -s ../kinetics_tools.tasks.ipd_summary-0/basemods.h5 output/tasks/kinetics_tools.tasks.gather_kinetics_h5-1/file.h5
