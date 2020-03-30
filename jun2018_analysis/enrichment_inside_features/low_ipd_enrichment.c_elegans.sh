#!/usr/bin/env bash
prefix=$(basename $0 .sh)
{
cat calc_enrichment_inside_features.csv_header
./calc_enrichment_inside_features.sh ../extreme_ipd.low.ab.coverage25.sorted.c_elegans.gff ../deep_kinetics_region.ab.bed.merged.sorted.c_elegans.bed ab
./calc_enrichment_inside_features.sh ../extreme_ipd.low.cd.coverage25.sorted.c_elegans.gff ../deep_kinetics_region.cd.bed.merged.sorted.c_elegans.bed cd
./calc_enrichment_inside_features.sh ../extreme_ipd.low.abcd.coverage25.sorted.c_elegans.gff ../deep_kinetics_region.abcd.bed.merged.sorted.c_elegans.bed abcd
./calc_enrichment_inside_features.sh ../extreme_ipd.low.PD2182.coverage25.sorted.c_elegans.gff ../deep_kinetics_region.PD2182.bed.merged.sorted.c_elegans.bed PD2182
./calc_enrichment_inside_features.sh ../extreme_ipd.low.PD2182sequel.coverage25.sorted.c_elegans.gff ../deep_kinetics_region.PD2182sequel.bed.merged.sorted.c_elegans.bed PD2182sequel
./calc_enrichment_inside_features.sh ../extreme_ipd.low.k.coverage25.sorted.c_elegans.gff ../deep_kinetics_region.k.bed.merged.sorted.c_elegans.bed k
./calc_enrichment_inside_features.sh ../extreme_ipd.low.l.coverage25.sorted.c_elegans.gff ../deep_kinetics_region.l.bed.merged.sorted.c_elegans.bed l
} > $prefix.csv 2> $prefix.log
