#!/usr/bin/env bash
MOTIF=$(cat PATTERN)
ecoli="E._coli"

output_name_celegans="C. elegans"
output_name_ecoli="E. coli"

# header
echo "motif,sample_id,feature,chromosome_set,num_occ"

ref_id="ce11rel606"
ref_feature_name="whole_genome"
ref_celegans=$(cat motif_occ.ce11rel606.bed.merged_occ | (grep -v "$ecoli" || [[ $? == 1 ]]) | wc -l)
ref_ecoli=$(cat motif_occ.ce11rel606.bed.merged_occ | (grep "$ecoli" || [[ $? == 1 ]]) | wc -l)
echo "$MOTIF,$ref_id,$ref_feature_name,$output_name_celegans,$ref_celegans"
echo "$MOTIF,$ref_id,$ref_feature_name,$output_name_ecoli,$ref_ecoli"

sample_ids="ab cd k l PD2182sequel abcd kl PD2182"
for ID in $sample_ids
do
    deep_feature_name="deep_kinetics_region"
    deep_celegans=$(cat motif_occ.deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.bed.merged_occ | (grep -v "$ecoli" || [[ $? == 1 ]]) | wc -l)
    deep_ecoli=$(cat motif_occ.deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.bed.merged_occ | (grep "$ecoli" || [[ $? == 1 ]]) | wc -l)
    echo "$MOTIF,$ID,$deep_feature_name,$output_name_celegans,$deep_celegans"
    echo "$MOTIF,$ID,$deep_feature_name,$output_name_ecoli,$deep_ecoli"

    high_feature_name="high IPD"
    high_celegans=$(cat motif_occ.extreme_ipd.high.$ID.coverage25.slop20.fullLength.gff.sorted.bed.bed.merged_occ | (grep -v "$ecoli" || [[ $? == 1 ]]) | wc -l)
    high_ecoli=$(cat motif_occ.extreme_ipd.high.$ID.coverage25.slop20.fullLength.gff.sorted.bed.bed.merged_occ | (grep "$ecoli" || [[ $? == 1 ]]) | wc -l)
    echo "$MOTIF,$ID,$high_feature_name,$output_name_celegans,$high_celegans"
    echo "$MOTIF,$ID,$high_feature_name,$output_name_ecoli,$high_ecoli"

    low_feature_name="low IPD"
    low_celegans=$(cat motif_occ.extreme_ipd.low.$ID.coverage25.slop20.fullLength.gff.sorted.bed.bed.merged_occ | (grep -v "$ecoli" || [[ $? == 1 ]]) | wc -l)
    low_ecoli=$(cat motif_occ.extreme_ipd.low.$ID.coverage25.slop20.fullLength.gff.sorted.bed.bed.merged_occ | (grep "$ecoli" || [[ $? == 1 ]]) | wc -l)
    echo "$MOTIF,$ID,$low_feature_name,$output_name_celegans,$low_celegans"
    echo "$MOTIF,$ID,$low_feature_name,$output_name_ecoli,$low_ecoli"
done
