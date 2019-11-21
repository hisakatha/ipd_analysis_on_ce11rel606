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

mod_csv=()
while read LINE
do
    mod_csv+=("$LINE")
done < <(tail -n +2 ../get_enriched_gene.both.summarize.csv)

for i in $(seq 0 $((NUM - 1)))
do
    SAMPLE1=${mod_csv[$((2 * i))]}
    SAMPLE2=${mod_csv[$((2 * i + 1))]}
    # csvcut is a part of csvkit, a python module
    GENE=$(echo $SAMPLE1 | csvcut -c 1)
    gene_class=$(curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/${gene_ids[$i]}/gene_class | jq .gene_class.data.description | sed -E 's/"(.*)"/\1/; t; s/^$/(NA)/')
    description=$(curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/${gene_ids[$i]}/concise_description | jq .concise_description.data.text | sed -E 's/"(.*)"/\1/; t; s/^$/(NA)/')
    FEATURE=$(echo $SAMPLE1 | csvcut -c 2 | sed -E 's/"(.*)"/\1/')
    callable1=$(echo $SAMPLE1 | csvcut -c 4)
    modified1=$(echo $SAMPLE1 | csvcut -c 5)
    fraction1=$(echo $SAMPLE1 | csvcut -c 6)
    pvalue1=$(echo $SAMPLE1 | csvcut -c 7)
    qvalue1=$(echo $SAMPLE1 | csvcut -c 10)
    callable2=$(echo $SAMPLE2 | csvcut -c 4)
    modified2=$(echo $SAMPLE2 | csvcut -c 5)
    fraction2=$(echo $SAMPLE2 | csvcut -c 6)
    pvalue2=$(echo $SAMPLE2 | csvcut -c 7)
    qvalue2=$(echo $SAMPLE2 | csvcut -c 10)
    echo "Gene: $GENE"; echo
    echo "Type: \verb|${biotypes[$i]}|"; echo
    echo "Common name: ${common_names[$i]}"; echo
    echo "Gene class: $gene_class"; echo
    echo "Description: $description"; echo
    echo '\medskip'
    echo "Enriched feature type: $FEATURE"; echo
    echo '\begin{tabular}{ l r r l l l }'
    echo 'Replicate & \#callable & \#modified & Fraction & p-value & FDR q-value \\'
    echo "Sample\\#1 & $callable1 & $modified1 & $fraction1 & $pvalue1 & $qvalue1"'\\'
    echo "Sample\\#2 & $callable2 & $modified2 & $fraction2 & $pvalue2 & $qvalue2"
    echo '\end{tabular}'
    echo '\medskip'
    echo '\begin{scriptsize}'
    echo '\color{gray}'
    #echo '\verb'"|${figures[$i]}|"; echo
    #echo '\verb'"|$(echo $SAMPLE1 | sed -E 's/,([^ ])/, \1/g')|"; echo
    #echo '\verb'"|$(echo $SAMPLE2 | sed -E 's/,([^ ])/, \1/g')|"; echo
    echo '\begin{lstlisting}[breaklines]'
    echo "${figures[$i]}"
    echo "$(echo $SAMPLE1 | sed -E 's/,([^ ])/, \1/g')"
    echo "$(echo $SAMPLE2 | sed -E 's/,([^ ])/, \1/g')"
    echo '\end{lstlisting}'
    echo '\end{scriptsize}'
    echo '\includegraphics[width=\linewidth]'"{{${figures[$i]}}}"
    echo '\newpage'
    echo "Processed $((i + 1)) / $NUM" >&2
done >> $OUT

echo '\end{document}' >> $OUT

lualatex $OUT
