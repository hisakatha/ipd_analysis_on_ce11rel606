#!/usr/bin/env bash
# This script tests the difference between the two sets in the following contingency table:
# (the number of occurrence of a motif in a feature region (high IPD or low IPD)) vs. (the number of occurrence of a motif in deep_kinetics_region), and
# (the number of bases in a feature region) vs. (the number of bases in deep_kinetics_region)
#
# deep_kinetics_region: region with a valid IPD count more than a threshold (default: count >= 25)

MOTIF=$(cat PATTERN)
# sample id (e.g. ab, cd)
ID=$1
ecoli="E._coli"
high_test_num=$2
low_test_num=$3

bg=$(cat $MOTIF.deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.fa.bed.merged_occ | (grep -v "$ecoli" || [[ $? == 1 ]]) | wc -l)
high=$(cat $MOTIF.extreme_ipd.high.$ID.coverage25.slop20.fullLength.gff.fa.bed.merged_occ | (grep -v "$ecoli" || [[ $? == 1 ]]) | wc -l)
low=$(cat $MOTIF.extreme_ipd.low.$ID.coverage25.slop20.fullLength.gff.fa.bed.merged_occ | (grep -v "$ecoli" || [[ $? == 1 ]]) | wc -l)
all_bg=$(cat ../../deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.merged | (grep -v "$ecoli" || [[ $? == 1 ]]) | awk '{sum += $3 - $2}END{print sum}')
all_high=$(cat ../../extreme_ipd.high.$ID.coverage25.slop20.fullLength.gff.merged | (grep -v "$ecoli" || [[ $? == 1 ]]) | awk '{sum += $3 - $2}END{print sum}')
all_low=$(cat ../../extreme_ipd.low.$ID.coverage25.slop20.fullLength.gff.merged | (grep -v "$ecoli" || [[ $? == 1 ]]) | awk '{sum += $3 - $2}END{print sum}')
bg_ecoli=$(cat $MOTIF.deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.fa.bed.merged_occ | (grep "$ecoli" || [[ $? == 1 ]]) | wc -l)
high_ecoli=$(cat $MOTIF.extreme_ipd.high.$ID.coverage25.slop20.fullLength.gff.fa.bed.merged_occ | (grep "$ecoli" || [[ $? == 1 ]]) | wc -l)
low_ecoli=$(cat $MOTIF.extreme_ipd.low.$ID.coverage25.slop20.fullLength.gff.fa.bed.merged_occ | (grep "$ecoli" || [[ $? == 1 ]]) | wc -l)
all_bg_ecoli=$(cat ../../deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.merged | (grep "$ecoli" || [[ $? == 1 ]]) | awk '{sum += $3 - $2}END{print sum}')
all_high_ecoli=$(cat ../../extreme_ipd.high.$ID.coverage25.slop20.fullLength.gff.merged | (grep "$ecoli" || [[ $? == 1 ]]) | awk '{sum += $3 - $2}END{print sum}')
all_low_ecoli=$(cat ../../extreme_ipd.low.$ID.coverage25.slop20.fullLength.gff.merged | (grep "$ecoli" || [[ $? == 1 ]]) | awk '{sum += $3 - $2}END{print sum}')
echo "# all_bg_celegans = $all_bg, all_high_celegans = $all_high, all_low_celegans = $all_low"
echo "# all_bg_ecoli = $all_bg_ecoli, all_high_ecoli = $all_high_ecoli, all_low_ecoli = $all_low_ecoli"

high_ratio=$(awk -v feature=$high -v bg=$bg 'BEGIN{printf "%.3f\n", feature / bg}')
high_p=$(Rscript -e "cat(sprintf('%.3g\n', chisq.test(matrix(c($all_high, $all_bg, $high, $bg), 2))\$p.value))")
high_adjusted_p=$(awk -v pvalue=$high_p -v test_num=$high_test_num 'BEGIN{printf "%.3g\n", pvalue * test_num}')
low_ratio=$(awk -v feature=$low -v bg=$bg 'BEGIN{printf "%.3f\n", feature / bg}')
low_p=$(Rscript -e "cat(sprintf('%.3g\n', chisq.test(matrix(c($all_low, $all_bg, $low, $bg), 2))\$p.value))")
low_adjusted_p=$(awk -v pvalue=$low_p -v test_num=$low_test_num 'BEGIN{printf "%.3g\n", pvalue * test_num}')
echo "sample_id,species,motif,feature_type,occ_feature,occ_background,ratio,chisq.pvalue,adjusted.pvalue"
echo "$ID,C. elegans,$MOTIF,high_ipd,$high,$bg,$high_ratio,$high_p,$high_adjusted_p"
echo "$ID,C. elegans,$MOTIF,low_ipd,$low,$bg,$low_ratio,$low_p,$low_adjusted_p"

high_ratio=$(awk -v feature=$high_ecoli -v bg=$bg_ecoli 'BEGIN{printf "%.3f\n", feature / bg}')
high_p=$(Rscript -e "cat(sprintf('%.3g\n', chisq.test(matrix(c($all_high_ecoli, $all_bg_ecoli, $high_ecoli, $bg_ecoli), 2))\$p.value))")
high_adjusted_p=$(awk -v pvalue=$high_p -v test_num=$high_test_num 'BEGIN{printf "%.3g\n", pvalue * test_num}')
low_ratio=$(awk -v feature=$low_ecoli -v bg=$bg_ecoli 'BEGIN{printf "%.3f\n", feature / bg}')
low_p=$(Rscript -e "cat(sprintf('%.3g\n', chisq.test(matrix(c($all_low_ecoli, $all_bg_ecoli, $low_ecoli, $bg_ecoli), 2))\$p.value))")
low_adjusted_p=$(awk -v pvalue=$low_p -v test_num=$low_test_num 'BEGIN{printf "%.3g\n", pvalue * test_num}')
#echo "sample_id,species,motif,feature_type,occ_feature,occ_background,ratio,chisq.pvalue,adjusted.pvalue"
echo "$ID,E. coli,$MOTIF,high_ipd,$high_ecoli,$bg_ecoli,$high_ratio,$high_p,$high_adjusted_p"
echo "$ID,E. coli,$MOTIF,low_ipd,$low_ecoli,$bg_ecoli,$low_ratio,$low_p,$low_adjusted_p"
