#!/usr/bin/env bash

write_header () {
    cat <<'_EOS_'
\documentclass[12pt,a4paper]{article}
\usepackage[margin=2mm]{geometry}
\usepackage{lmodern}
\usepackage[T1]{fontenc}
\usepackage[hiresbb]{graphicx}
\usepackage{grffile}
\begin{document}
_EOS_
}

dir_to_title () {
    if [[ $1 == "GGN_4" || $1 == "GGN_10" ]]
    then
        title=$(echo $1 | sed -E "s/([a-zA-Z]+)_([0-9]+)/(\1)\2/")
    else
        title=$1
    fi
    echo $title
}

write_include_for_motif () {
    original_pdf="$1/plot_kinetics_correlation.celegans_ecoli.pdf"
    cropped_pdf="cropped_kinetics_correlation.celegans_ecoli.$1.pdf"
    pdfcrop --hires "$original_pdf" "$cropped_pdf" >&2
    echo '\begin{minipage}[t]{0.33\linewidth}'
    title=$(dir_to_title $1)
    echo '\verb'"|${title}|"'\\'
    echo '\includegraphics[page=2,width=0.95\linewidth]'"{{$cropped_pdf}}"
    echo '\end{minipage}'
}

write_footer () {
    echo '\end{document}'
}

write_all_for_motifs () {
    # Accept N arguments (N >= 2)
    # 1st: output file name
    # 2nd to Nth: (N - 1) motif directory names
    OUT="$1"
    shift
    MOTIFS=("$@")
    for i in "${MOTIFS[@]}"
    do
        write_include_for_motif "$i"
    done >> $OUT
}

#motifs_high_ipd="ACGCRTG ATCAGCTG GGN_4 RGTA AAAABT AGAGAGTA AGCTATAT AGGCAGGC ATGGGAYA ATTGTTAC CAAACTAC CAGYTG CCAATCAG CRACGAS DCGAGACC DGCTTC GAAGGATC GATATRGY GCGACCTA GCGCGCGC GGHGGY GTAGATCA GTATCGTA TGACGTCA TGGTGSA YGGAR"
motifs_high_ipd1="ACGCRTG ATCAGCTG GGN_4 RGTA AAAABT AGAGAGTA AGCTATAT AGGCAGGC ATGGGAYA"
motifs_high_ipd2="ATTGTTAC CAAACTAC CAGYTG CCAATCAG CRACGAS DCGAGACC DGCTTC GAAGGATC GATATRGY"
motifs_high_ipd3="GCGACCTA GCGCGCGC GGHGGY GTAGATCA GTATCGTA TGACGTCA TGGTGSA YGGAR"

motifs_low_ipd1="ACATMTGG CTGDAR TACTGTAG AATMAATA AGACGCAG CGGYTTGA CTKCAA GCGCGTCA GTGTGTGY"
motifs_low_ipd2="TACCCCKA TACCTTGA TGACGTCA"

TEX="arrange_plot_kinetics_correlation.extreme_ipd.celegans_ecoli_v2.tex"
echo "" > $TEX
write_header >> $TEX
# High IPD motifs
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs_high_ipd1
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs_high_ipd2
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs_high_ipd3
echo '\end{figure}' >> $TEX
# Low IPD motifs
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs_low_ipd1
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs_low_ipd2
echo '\end{figure}' >> $TEX
write_footer >> $TEX

lualatex $TEX
pdfcrop --hires "$(basename $TEX .tex).pdf"
