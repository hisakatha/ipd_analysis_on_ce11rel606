#!/usr/bin/env bash
if [[ $# < 2 || $# > 3 ]]; then
    echo "Usage: $0 GENOME INCLUDE_BED INPUT" >&2
    echo "    GENOME: tab-deliminated chromInfo or FASTA index" >&2
    echo "    INCLUDE_BED: regions in which shuffled start position should be placed in BED format" >&2
    echo "    INPUT: regions to be shuffled in BED/GFF format. GFF files need to have .gff[3] extension, otherwise input will be regarded as BED." >&2
    exit 1
fi

genome="$1"
incl="$2"
input="$3"

# temp files for using starnd information
# TODO?: assign a fixed name if process substitution is used.
# Maybe, fixed names are bad for parallel invocation.
#if [[ -p $1 ]]; then echo $1 is a pipe && genome_s=GENOME.with_strand; fi; exit
#genome_s=$(basename $genome).with_strand
#incl_s=$(basename $incl).with_strand
genome_s=$(mktemp $(basename $genome).with_strand.XXXXXX)
incl_s=$(mktemp $(basename $incl).with_strand.XXXXXX)

cat $genome | sed -E 's/^([^ \t]+)/\1_pos/' > $genome_s
cat $genome | sed -E 's/^([^ \t]+)/\1_neg/' >> $genome_s

# TODO: handle an undefined strand "."
# TODO: handle GFF/VCF incl
cat $incl | awk -v FS="\t" -v OFS="\t" '$6 == "+"{print $1 "_pos",$2,$3} $6 == "-"{print $1 "_neg",$2,$3}' > $incl_s

# TODO: handle VCF input
SEED="-seed 1234"
fileType=$(echo $input | sed -E "s/^.*\.gff(3)?$/GFF/; t; s/.*/BED/")
cat $input | bedtools shuffle -g $genome_s -incl $incl_s $SEED |
    awk -v FS="\t" -v OFS="\t" -v fileType=$fileType \
    '{if(match($1, "_pos$") > 0){$1 = substr($1, 1, length($1) - 4); strand = "+"} else if(match($1, "_neg$") > 0){$1 = substr($1, 1, length($1) - 4); strand = "-"} else {strand = "."}; if(fileType=="GFF"){$7 = strand; print $0}else{print $1,$2,$3,".",".",strand}}'

rm $genome_s $incl_s

