#!/usr/bin/env bash
MOTIF=$(cat PATTERN)
# sample id (e.g. ab, cd)
ID=$1
ecoli="E._coli"
high_test_num=8
low_test_num=6

bg=$(cat $MOTIF.deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.fa.bed.merged_occ | (grep -v "$ecoli" || [[ $? == 1 ]]) | wc -l)
high=$(cat $MOTIF.extreme_ipd.high.$ID.coverage25.slop20.fullLength.gff.fa.bed.merged_occ | (grep -v "$ecoli" || [[ $? == 1 ]]) | wc -l)
low=$(cat $MOTIF.extreme_ipd.low.$ID.coverage25.slop20.fullLength.gff.fa.bed.merged_occ | (grep -v "$ecoli" || [[ $? == 1 ]]) | wc -l)
all_bg=$(cat ../../deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.merged | (grep -v "$ecoli" || [[ $? == 1 ]]) | awk '{sum += $3 - $2}END{print sum}')
all_high=$(cat ../../extreme_ipd.high.$ID.coverage25.slop20.fullLength.gff.merged | (grep -v "$ecoli" || [[ $? == 1 ]]) | awk '{sum += $3 - $2}END{print sum}')
all_low=$(cat ../../extreme_ipd.low.$ID.coverage25.slop20.fullLength.gff.merged | (grep -v "$ecoli" || [[ $? == 1 ]]) | awk '{sum += $3 - $2}END{print sum}')
echo "# all_bg = $all_bg, all_high = $all_high, all_low = $all_low"
high_ratio=$(awk -v feature=$high -v bg=$bg 'BEGIN{printf "%.3f\n", feature / bg}')
high_p=$(Rscript -e "cat(sprintf('%.3g\n', chisq.test(matrix(c($all_high, $all_bg, $high, $bg), 2))\$p.value))")
high_adjusted_p=$(awk -v pvalue=$high_p -v test_num=$high_test_num 'BEGIN{printf "%.3g\n", pvalue * test_num}')
low_ratio=$(awk -v feature=$low -v bg=$bg 'BEGIN{printf "%.3f\n", feature / bg}')
low_p=$(Rscript -e "cat(sprintf('%.3g\n', chisq.test(matrix(c($all_low, $all_bg, $low, $bg), 2))\$p.value))")
low_adjusted_p=$(awk -v pvalue=$low_p -v test_num=$low_test_num 'BEGIN{printf "%.3g\n", pvalue * test_num}')
echo "sample_id,motif,feature_type,occ_feature,occ_background,ratio,chisq.pvalue,adjusted.pvalue"
echo "$ID,$MOTIF,high_ipd,$high,$bg,$high_ratio,$high_p,$high_adjusted_p"
echo "$ID,$MOTIF,low_ipd,$low,$bg,$low_ratio,$low_p,$low_adjusted_p"
