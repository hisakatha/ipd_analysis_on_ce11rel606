#!/usr/bin/env bash
prefix=$(basename $0 .sh)
{
cat calc_enrichment_inside_features.csv_header
./calc_enrichment_inside_features.sh ../extreme_ipd.low.ab.coverage25.sorted.gff ../deep_kinetics_region.ab.bed.merged.sorted ab
./calc_enrichment_inside_features.sh ../extreme_ipd.low.cd.coverage25.sorted.gff ../deep_kinetics_region.cd.bed.merged.sorted cd
./calc_enrichment_inside_features.sh ../extreme_ipd.low.abcd.coverage25.sorted.gff ../deep_kinetics_region.abcd.bed.merged.sorted abcd
./calc_enrichment_inside_features.sh ../extreme_ipd.low.PD2182.coverage25.sorted.gff ../deep_kinetics_region.PD2182.bed.merged.sorted PD2182
./calc_enrichment_inside_features.sh ../extreme_ipd.low.PD2182sequel.coverage25.sorted.gff ../deep_kinetics_region.PD2182sequel.bed.merged.sorted PD2182sequel
} > $prefix.csv 2> $prefix.log
