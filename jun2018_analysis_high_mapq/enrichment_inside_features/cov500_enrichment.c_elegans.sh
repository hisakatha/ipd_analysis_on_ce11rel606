#!/usr/bin/env bash
prefix=$(basename $0 .sh)
{
cat calc_enrichment_inside_features.csv_header
./calc_enrichment_inside_features.bed_input.sh ../deep_kinetics_region_cov500.ab.bed.merged.sorted.c_elegans.bed ../deep_kinetics_region.ab.bed.merged.sorted.c_elegans.bed ab
./calc_enrichment_inside_features.bed_input.sh ../deep_kinetics_region_cov500.cd.bed.merged.sorted.c_elegans.bed ../deep_kinetics_region.cd.bed.merged.sorted.c_elegans.bed cd
./calc_enrichment_inside_features.bed_input.sh ../deep_kinetics_region_cov500.PD2182.bed.merged.sorted.c_elegans.bed ../deep_kinetics_region.PD2182.bed.merged.sorted.c_elegans.bed PD2182
./calc_enrichment_inside_features.bed_input.sh ../deep_kinetics_region_cov500.PD2182sequel.bed.merged.sorted.c_elegans.bed ../deep_kinetics_region.PD2182sequel.bed.merged.sorted.c_elegans.bed PD2182sequel
./calc_enrichment_inside_features.bed_input.sh ../deep_kinetics_region_cov500.k.bed.merged.sorted.c_elegans.bed ../deep_kinetics_region.k.bed.merged.sorted.c_elegans.bed k
./calc_enrichment_inside_features.bed_input.sh ../deep_kinetics_region_cov500.l.bed.merged.sorted.c_elegans.bed ../deep_kinetics_region.l.bed.merged.sorted.c_elegans.bed l
} > $prefix.csv 2> $prefix.log
