#!/bin/bash

# INPUT: PacBio Sequel *.subreads.bam files
/bio/package/pacbio/smrtlink/smrtcmds/bin/dataset create --type SubreadSet "$(basename $(pwd)).subreads.xml" "$@"

