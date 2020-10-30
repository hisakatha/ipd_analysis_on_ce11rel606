#!/usr/bin/env bash
# Input regions in BED format
INPUT="$1"
INPUT_BASE="$(basename $INPUT .bed)"
# Callable region with deep coverage depths (sorted BED file)
# e.g. mapped.alignmentset.merged.bam.cov25
DEEP="$2"
FA="/glusterfs/hisakatha/ce11rel606/ce11rel606.fa"
FAI="/glusterfs/hisakatha/ce11rel606/ce11rel606.fa.fai"
# Input sample name without spaces or commas
INPUT_ID="$3"

get_bed_base_nums () {
    # Get the numbers of A, C, G, T and all the bases in a BED file
    bedtools nuc -s -fi $FA -bed $1 | grep -v "^#" | awk -v FS="\t" 'BEGIN{a=0;c=0;g=0;t=0} {a+=$9; c+=$10; g+=$11; t+=$12} END{print a,c,g,t,a+c+g+t}'
}

#get_gff_base_nums () {
#    # Get the numbers of A, C, G, T and all the bases in a GFF file
#    bedtools nuc -s -fi $FA -bed $1 | grep -v "^#" | awk -v FS="\t" 'BEGIN{a=0;c=0;g=0;t=0} {a+=$12; c+=$13; g+=$14; t+=$15} END{print a,c,g,t,a+c+g+t}'
#}

all_input_nums=($(get_bed_base_nums $INPUT))
all_deep_nums=($(get_bed_base_nums $DEEP))
base_names=("A" "C" "G" "T" "ALL")
all_fractions=()
for i in $(seq 0 4); do
    all_fractions+=($(awk -v denom=${all_deep_nums[$i]} -v numer=${all_input_nums[$i]} 'BEGIN{if (denom == 0) {print -1} else {printf "%.7g", numer / denom}}'))
done

# Features
FEATURE_BASE="/glusterfs/hisakatha/wormbase/c_elegans.PRJNA13758.WS267.annotations.gff3.gz.ucsc_chr"
ENHANCER="$FEATURE_BASE.enhancer.sorted.merged"
EXON="$FEATURE_BASE.exon.sorted.merged"
FIVE_PRIME_UTR="$FEATURE_BASE.five_prime_utr.sorted.merged"
INTRON="$FEATURE_BASE.intron.sorted.merged"
NCRNA="$FEATURE_BASE.ncrna.sorted.merged"
PROMOTER="$FEATURE_BASE.promoter.sorted.merged"
PSEUDOGENE_EXON="$FEATURE_BASE.pseudogene_exon.sorted.merged"
PSEUDOGENE="$FEATURE_BASE.pseudogene.sorted.merged"
TANDEM_REPEAT="$FEATURE_BASE.tandem_repeat.sorted.merged"
THREE_PRIME_UTR="$FEATURE_BASE.three_prime_utr.sorted.merged"
TSS_UP500="$FEATURE_BASE.tss_up500.sorted.merged"
TSS_DOWN500="$FEATURE_BASE.tss_down500.sorted.merged"
TTS_UP500="$FEATURE_BASE.tts_up500.sorted.merged"
TTS_DOWN500="$FEATURE_BASE.tts_down500.sorted.merged"
REPEAT_REGION="$FEATURE_BASE.repeat_region.sorted.merged"

calc_enrichment_in_a_feature () {
    FEATURE_FILE=$1
    FEATURE_ID=$2
    FEATURE_NAME="$3"
    echo "Start calculating enrichment in $FEATURE_ID" >&2
    # "-s" (force the same strands) or "" (no strands)
    STRAND_FLAG=""
    # You can skip sorting if all the input regions have 1 bp length
    bedtools intersect $STRAND_FLAG -sorted -g $FAI -a $INPUT -b $FEATURE_FILE | bedtools sort -faidx $FAI | bedtools merge -s -c 4,5,6 -o distinct,mean,distinct > $INPUT_BASE.inside.$FEATURE_ID
    bedtools intersect $STRAND_FLAG -sorted -g $FAI -a $DEEP -b $FEATURE_FILE | bedtools sort -faidx $FAI | bedtools merge -s -c 4,5,6 -o distinct,mean,distinct > deep_region.$INPUT_ID.inside.$FEATURE_ID
    #input_with_feature=$(cat $INPUT_BASE.inside.$FEATURE_ID | wc -l)
    input_with_feature_nums=($(get_bed_base_nums $INPUT_BASE.inside.$FEATURE_ID))
    #callable_with_feature=$(cat deep_region.$INPUT_ID.inside.$FEATURE_ID | awk -v total=0 '{total += $3 - $2} END{print total}')
    callable_with_feature_nums=($(get_bed_base_nums deep_region.$INPUT_ID.inside.$FEATURE_ID))
    for i in $(seq 0 4); do
        fraction=$(awk -v denom=${callable_with_feature_nums[$i]} -v numer=${input_with_feature_nums[$i]} 'BEGIN{if (denom == 0) {print -1} else {printf "%.7g", numer / denom}}')
        enrichment=$(awk -v denom=${all_fractions[$i]} -v numer=${fraction} 'BEGIN{if (denom == 0) {print -1} else {printf "%.7g", numer / denom}}')
        pvalue=$(Rscript -e "cat(binom.test(${input_with_feature_nums[$i]}, ${callable_with_feature_nums[$i]}, p = ${all_fractions[$i]})\$p.value)")
        echo "$INPUT_ID,\"$FEATURE_NAME\",${base_names[$i]},${callable_with_feature_nums[$i]},${input_with_feature_nums[$i]},$fraction,$enrichment,$pvalue"
    done
}

#echo "sample,region,base,callableNum,modifiedNum,fraction,enrichment,pvalue"
for i in $(seq 0 4); do
    echo "$INPUT_ID,ALL,${base_names[$i]},${all_deep_nums[$i]},${all_input_nums[$i]},${all_fractions[$i]},1,1"
done

calc_enrichment_in_a_feature $ENHANCER enhancer enhancer
calc_enrichment_in_a_feature $EXON exon exon
calc_enrichment_in_a_feature $FIVE_PRIME_UTR five_prime_utr "five prime utr"
calc_enrichment_in_a_feature $INTRON intron intron
calc_enrichment_in_a_feature $NCRNA ncRNA ncRNA
calc_enrichment_in_a_feature $PROMOTER promoter promoter
calc_enrichment_in_a_feature $PSEUDOGENE_EXON pseudogene_exon "pseudogene exon"
calc_enrichment_in_a_feature $PSEUDOGENE pseudogene "pseudogene"
calc_enrichment_in_a_feature $TANDEM_REPEAT tandem_repeat "tandem repeat"
calc_enrichment_in_a_feature $THREE_PRIME_UTR three_prime_utr "three prime utr"
calc_enrichment_in_a_feature $TSS_UP500 tss_up500 "TSS [-500, 0]"
calc_enrichment_in_a_feature $TSS_DOWN500 tss_down500 "TSS [0, 500]"
calc_enrichment_in_a_feature $TTS_UP500 tts_up500 "TTS [-500, 0]"
calc_enrichment_in_a_feature $TTS_DOWN500 tts_down500 "TTS [0, 500]"
calc_enrichment_in_a_feature $REPEAT_REGION repeat_region "repeat region"

