#!/usr/bin/env bash
# TODO?: Fisher's exact test assumes fixed merginal counts, which is not the case. We need to find another test?
# This script tests the difference between the two sets in the following contingency table:
# (the number of occurrence of a motif in a feature region (high IPD or low IPD)) vs. (the number of occurrence of a motif outside the feature region), and
# (the number of bases in a feature region) vs. (the number of bases outside the motif region)
#
# deep_kinetics_region: region with a valid IPD count more than a threshold (default: count >= 25)

set -u
MOTIF=$(cat PATTERN)
# sample id (e.g. ab, cd)
ID=$1
ecoli="E._coli"
high_test_num=$2
low_test_num=$3
output_name_celegans="C. elegans"
output_name_ecoli="E. coli"

bg=$(cat $MOTIF.deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.fa.bed.merged_occ | (grep -v "$ecoli" || [[ $? == 1 ]]) | wc -l)
high=$(cat $MOTIF.extreme_ipd.high.$ID.coverage25.slop20.fullLength.gff.fa.bed.merged_occ | (grep -v "$ecoli" || [[ $? == 1 ]]) | wc -l)
nonhigh=$((bg - high))
low=$(cat $MOTIF.extreme_ipd.low.$ID.coverage25.slop20.fullLength.gff.fa.bed.merged_occ | (grep -v "$ecoli" || [[ $? == 1 ]]) | wc -l)
nonlow=$((bg - low))
all_bg=$(cat ../../deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.merged | (grep -v "$ecoli" || [[ $? == 1 ]]) | awk 'BEGIN{sum = 0}{sum += $3 - $2}END{print sum}')
all_high=$(cat ../../extreme_ipd.high.$ID.coverage25.slop20.fullLength.gff.merged | (grep -v "$ecoli" || [[ $? == 1 ]]) | awk 'BEGIN{sum = 0}{sum += $3 - $2}END{print sum}')
all_nonhigh=$((all_bg - all_high))
all_low=$(cat ../../extreme_ipd.low.$ID.coverage25.slop20.fullLength.gff.merged | (grep -v "$ecoli" || [[ $? == 1 ]]) | awk 'BEGIN{sum = 0}{sum += $3 - $2}END{print sum}')
all_nonlow=$((all_bg - all_low))
bg_ecoli=$(cat $MOTIF.deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.fa.bed.merged_occ | (grep "$ecoli" || [[ $? == 1 ]]) | wc -l)
high_ecoli=$(cat $MOTIF.extreme_ipd.high.$ID.coverage25.slop20.fullLength.gff.fa.bed.merged_occ | (grep "$ecoli" || [[ $? == 1 ]]) | wc -l)
nonhigh_ecoli=$((bg_ecoli - high_ecoli))
low_ecoli=$(cat $MOTIF.extreme_ipd.low.$ID.coverage25.slop20.fullLength.gff.fa.bed.merged_occ | (grep "$ecoli" || [[ $? == 1 ]]) | wc -l)
nonlow_ecoli=$((bg_ecoli - low_ecoli))
all_bg_ecoli=$(cat ../../deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.merged | (grep "$ecoli" || [[ $? == 1 ]]) | awk 'BEGIN{sum = 0}{sum += $3 - $2}END{print sum}')
all_high_ecoli=$(cat ../../extreme_ipd.high.$ID.coverage25.slop20.fullLength.gff.merged | (grep "$ecoli" || [[ $? == 1 ]]) | awk 'BEGIN{sum = 0}{sum += $3 - $2}END{print sum}')
all_nonhigh_ecoli=$((all_bg_ecoli - all_high_ecoli))
all_low_ecoli=$(cat ../../extreme_ipd.low.$ID.coverage25.slop20.fullLength.gff.merged | (grep "$ecoli" || [[ $? == 1 ]]) | awk 'BEGIN{sum = 0}{sum += $3 - $2}END{print sum}')
all_nonlow_ecoli=$((all_bg_ecoli - all_low_ecoli))
echo "# all_bg_celegans = $all_bg, all_high_celegans = $all_high, all_low_celegans = $all_low"
echo "# all_bg_ecoli = $all_bg_ecoli, all_high_ecoli = $all_high_ecoli, all_low_ecoli = $all_low_ecoli"

high_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $high / $bg))")
all_high_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $all_high / $all_bg))")
high_odds_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $high / $nonhigh / $all_high * $all_nonhigh))")
high_p=$(python3 -c "import scipy.stats as ss; oddsratio,pvalue = ss.fisher_exact([[$all_high, $all_nonhigh], [$high, $nonhigh]]); print(pvalue)")
high_adjusted_p=$(Rscript -e "cat(sprintf('%.3g\n', $high_p * $high_test_num))")
echo "sample_id,chr,motif,feature_type,occ_feature,occ_nonfeature,occ_ratio,base_feature,base_nonfeature,base_ratio,odds_ratio,fisher.pvalue,adjusted.pvalue"
echo "$ID,$output_name_celegans,$MOTIF,high_ipd,$high,$nonhigh,$high_ratio,$all_high,$all_nonhigh,$all_high_ratio,$high_odds_ratio,$high_p,$high_adjusted_p"

low_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $low / $bg))")
all_low_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $all_low / $all_bg))")
low_odds_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $low / $nonlow / $all_low * $all_nonlow))")
low_p=$(python3 -c "import scipy.stats as ss; oddsratio,pvalue = ss.fisher_exact([[$all_low, $all_nonlow], [$low, $nonlow]]); print(pvalue)")
low_adjusted_p=$(Rscript -e "cat(sprintf('%.3g\n', $low_p * $low_test_num))")
echo "$ID,$output_name_celegans,$MOTIF,low_ipd,$low,$nonlow,$low_ratio,$all_low,$all_nonlow,$all_low_ratio,$low_odds_ratio,$low_p,$low_adjusted_p"

high_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $high_ecoli / $bg_ecoli))")
all_high_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $all_high_ecoli / $all_bg_ecoli))")
high_odds_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $high_ecoli / $nonhigh_ecoli / $all_high_ecoli * $all_nonhigh_ecoli))")
high_p=$(python3 -c "import scipy.stats as ss; oddsratio,pvalue = ss.fisher_exact([[$all_high_ecoli, $all_nonhigh_ecoli], [$high_ecoli, $nonhigh_ecoli]]); print(pvalue)")
high_adjusted_p=$(Rscript -e "cat(sprintf('%.3g\n', $high_p * $high_test_num))")
echo "$ID,$output_name_ecoli,$MOTIF,high_ipd,$high_ecoli,$nonhigh_ecoli,$high_ratio,$all_high_ecoli,$all_nonhigh_ecoli,$all_high_ratio,$high_odds_ratio,$high_p,$high_adjusted_p"

low_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $low_ecoli / $bg_ecoli))")
all_low_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $all_low_ecoli / $all_bg_ecoli))")
low_odds_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $low_ecoli / $nonlow_ecoli / $all_low_ecoli * $all_nonlow_ecoli))")
low_p=$(python3 -c "import scipy.stats as ss; oddsratio,pvalue = ss.fisher_exact([[$all_low_ecoli, $all_nonlow_ecoli], [$low_ecoli, $nonlow_ecoli]]); print(pvalue)")
low_adjusted_p=$(Rscript -e "cat(sprintf('%.3g\n', $low_p * $low_test_num))")
echo "$ID,$output_name_ecoli,$MOTIF,low_ipd,$low_ecoli,$nonlow_ecoli,$low_ratio,$all_low_ecoli,$all_nonlow_ecoli,$all_low_ratio,$low_odds_ratio,$low_p,$low_adjusted_p"
