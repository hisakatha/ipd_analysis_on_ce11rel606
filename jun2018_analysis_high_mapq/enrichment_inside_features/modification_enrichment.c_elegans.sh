#!/usr/bin/env bash
prefix=$(basename $0 .sh)
{
cat calc_enrichment_inside_features.csv_header
./calc_enrichment_inside_features.sh ../modifications.k_normBy_ab.score20_ipdratio2_coverage25.c_elegans.gff ../deep_kinetics_region.k_normBy_ab.bed.c_elegans.bed k_normBy_ab
./calc_enrichment_inside_features.sh ../modifications.l_normBy_cd.score20_ipdratio2_coverage25.c_elegans.gff ../deep_kinetics_region.l_normBy_cd.bed.c_elegans.bed l_normBy_cd
#./calc_enrichment_inside_features.sh ../modifications.kl_normBy_abcd.score20_ipdratio2_coverage25.c_elegans.gff ../deep_kinetics_region.kl_normBy_abcd.bed.c_elegans.bed kl_normBy_abcd
} > $prefix.csv 2> $prefix.log
