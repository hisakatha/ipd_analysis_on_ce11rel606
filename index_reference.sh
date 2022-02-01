#!/usr/bin/env bash
ref="/glusterfs/hisakatha/ce11rel606/ce11rel606.fa"
/bio/package/pacbio/smrtlink10/smrtcmds/bin/pbmm2 index --preset CCS "$ref" ce11rel606.ccs.mmi
/bio/package/pacbio/smrtlink10/smrtcmds/bin/pbmm2 index --preset SUBREAD "$ref" ce11rel606.subread.mmi
