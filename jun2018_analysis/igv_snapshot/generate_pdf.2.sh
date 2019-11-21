#!/usr/bin/env bash
# #!/bin/bash

OUT=get_enriched_genes.both.igv_snapshot.tex
cat <<'_EOS_' > $OUT
\documentclass[11pt,a4paper]{article}
\usepackage[margin=10mm,footskip=0.25in]{geometry}
\usepackage{lmodern}
\usepackage[T1]{fontenc}
\usepackage[hiresbb]{graphicx}
\usepackage{grffile}
\usepackage{listings}
\usepackage{xcolor}
\begin{document}
_EOS_

figures=($(cat ../get_enriched_gene.both.summarize.gene_region.igv_script | grep .svg | sed -e 's/snapshot //' -e 's/.svg/.png/'))

common_names=($(cat ../get_enriched_gene.both.summarize.gene_region.gff | cut -f9 | sed -E 's/.*locus=([-\._a-zA-Z0-9]+).*/\1/; t; s/.*/(NA)/'))
biotypes=($(cat ../get_enriched_gene.both.summarize.gene_region.gff | cut -f9 | sed -E 's/.*biotype=([-\._a-zA-Z0-9]+).*/\1/; t; s/.*/(NA)/'))
gene_ids=($(cat ../get_enriched_gene.both.summarize.gene_region.gff | cut -f9 | sed -E 's/.*ID=Gene:([-\._a-zA-Z0-9]+).*/\1/; t; s/.*/(NA)/'))

NUM=${#figures[*]}
#NUM=3

# csvcut is a part of csvkit, a python module
features=()
while read LINE
do
    features+=("$LINE")
done < <(tail -n +2 ../get_enriched_gene.both.summarize.csv | csvcut -c 2 | sed -E 's/"(.*)"/\1/')

genes=($(tail -n +2 ../get_enriched_gene.both.summarize.csv | csvcut -c 1))
samples=($(tail -n +2 ../get_enriched_gene.both.summarize.csv | csvcut -c 3 | sed -E 's/k_normBy_ab/VC2010+OP50/; s/l_normBy_cd/VC2010/'))
callable_nums=($(tail -n +2 ../get_enriched_gene.both.summarize.csv | csvcut -c 4))
modified_nums=($(tail -n +2 ../get_enriched_gene.both.summarize.csv | csvcut -c 5))
fraction_vals=($(tail -n +2 ../get_enriched_gene.both.summarize.csv | csvcut -c 6 | xargs -n1 printf "%.3g\n"))
pvalues=($(tail -n +2 ../get_enriched_gene.both.summarize.csv | csvcut -c 7 | xargs -n1 printf "%.3g\n"))
qvalues=($(tail -n +2 ../get_enriched_gene.both.summarize.csv | csvcut -c 10 | xargs -n1 printf "%.3g\n"))

sample1_all_frac=$(cat ../calc_enrichment.both.csv | csvgrep -c 1 -r "^k_normBy_ab$" | csvgrep -c 2 -r "^All$" | csvcut -c 5 | tail -n1)
sample2_all_frac=$(cat ../calc_enrichment.both.csv | csvgrep -c 1 -r "^l_normBy_cd$" | csvgrep -c 2 -r "^All$" | csvcut -c 5 | tail -n1)

for i in $(seq 0 $((NUM - 1)))
do
    i1=$((2 * i))
    i2=$((2 * i + 1))
    GENE=${genes[$i1]}
    gene_class=$(curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/${gene_ids[$i]}/gene_class | jq .gene_class.data.description | sed -E 's/"(.*)"/\1/; t; s/^$/(NA)/')
    description=$(curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/${gene_ids[$i]}/concise_description | jq .concise_description.data.text | sed -E 's/"(.*)"/\1/; t; s/^$/(NA)/')
    FEATURE=${features[$i1]}
    sample1=${samples[$i1]}
    callable1=${callable_nums[$i1]}
    modified1=${modified_nums[$i1]}
    fraction1=${fraction_vals[$i1]}
    enrichment1=$(echo "scale=4; $fraction1 / $sample1_all_frac" | bc | xargs -n1 printf "%.3g")
    pvalue1=${pvalues[$i1]}
    qvalue1=${qvalues[$i1]}
    sample2=${samples[$i2]}
    callable2=${callable_nums[$i2]}
    modified2=${modified_nums[$i2]}
    fraction2=${fraction_vals[$i2]}
    enrichment2=$(echo "scale=4; $fraction2 / $sample2_all_frac" | bc | xargs -n1 printf "%.3g")
    pvalue2=${pvalues[$i2]}
    qvalue2=${qvalues[$i2]}
    echo "Gene: $GENE"; echo
    echo "Type: \verb|${biotypes[$i]}|"; echo
    echo "Common name: ${common_names[$i]}"; echo
    echo "Gene class: $gene_class"; echo
    echo "Description: $description"; echo
    echo '\medskip'
    echo "Enriched feature type: $FEATURE"; echo
    echo '\begin{tabular}{ l r r l l l l }'
    echo 'Sample & \#callable & \#modified & Fraction & Enrichment & p-value & FDR q-value \\'
    echo "$sample1 & $callable1 & $modified1 & $fraction1 & $enrichment1 & $pvalue1 & $qvalue1"'\\'
    echo "$sample2 & $callable2 & $modified2 & $fraction2 & $enrichment2 & $pvalue2 & $qvalue2"
    echo '\end{tabular}'
    echo '\medskip'
    echo '\begin{scriptsize}'
    echo '\color{gray}'
    #echo '\verb'"|${figures[$i]}|"; echo
    #echo '\verb'"|$(echo $SAMPLE1 | sed -E 's/,([^ ])/, \1/g')|"; echo
    #echo '\verb'"|$(echo $SAMPLE2 | sed -E 's/,([^ ])/, \1/g')|"; echo
    echo '\begin{lstlisting}[breaklines]'
    echo "${figures[$i]}"
    #echo "$(echo $SAMPLE1 | sed -E 's/,([^ ])/, \1/g')"
    #echo "$(echo $SAMPLE2 | sed -E 's/,([^ ])/, \1/g')"
    echo '\end{lstlisting}'
    echo '\end{scriptsize}'
    echo '\includegraphics[width=\linewidth]'"{{${figures[$i]}}}"
    echo '\newpage'
    echo "Processed $((i + 1)) / $NUM" >&2
done >> $OUT

echo '\end{document}' >> $OUT

lualatex $OUT
