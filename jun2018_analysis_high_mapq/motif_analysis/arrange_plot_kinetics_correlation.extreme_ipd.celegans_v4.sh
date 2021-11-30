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
    original_pdf="$1/plot_kinetics_correlation.celegans_ecoli_subset2_positive_strand.pdf"
    cropped_pdf="cropped_kinetics_correlation.celegans_ecoli_subset2_positive_strand.$1.pdf"
    pdfcrop --hires "$original_pdf" "$cropped_pdf" >&2
    echo '\begin{minipage}[t]{0.33\linewidth}'
    title=$(dir_to_title $1)
    echo '\verb'"|${title}|"'\\'
    echo '\includegraphics[page=8,width=0.95\linewidth]'"{{$cropped_pdf}}"
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

#motifs_high_ipd1="ACGCRTG ATCAGCTG GGN_4 RGTA AAAABT AGAGAGTA AGCTATAT AGGCAGGC ATGGGAYA ATTGTTAC CAAACTAC CAGYTG CCAATCAG CRACGAS DCGAGACC DGCTTC"
#motifs_high_ipd2="GAAGGATC GATATRGY GCGACCTA GCGCGCGC GGHGGY GTAGATCA GTATCGTA TGACGTCA TGGTGSA YGGAR"
#motifs_low_ipd="ACATMTGG CTGDAR TACTGTAG AATMAATA AGACGCAG CGGYTTGA CTKCAA GCGCGTCA GTGTGTGY TACCCCKA TACCTTGA TGACGTCA"

motifs1="ACGCRTG ATCAGCTG GGN_4"
motifs2="RGTA AAAABT AGAGAGTA"
motifs3="AGCTATAT AGGCAGGC ATGGGAYA"
motifs4="ATTGTTAC CAAACTAC CAGYTG"
motifs5="CCAATCAG CRACGAS DCGAGACC"
motifs6="DGCTTC GAAGGATC GATATRGY"
motifs7="GCGACCTA GCGCGCGC GGHGGY"
motifs8="GTAGATCA GTATCGTA TGACGTCA"
#motifs9="TGGTGSA YGGAR ACATMTGG"
#motifs10="CTGDAR TACTGTAG AATMAATA"
#motifs11="AGACGCAG CGGYTTGA CTKCAA"
#motifs12="GCGCGTCA GTGTGTGY TACCCCKA"
#motifs13="TACCTTGA"
motifs9="TGGTGSA YGGAR"
motifs10="ACATMTGG CTGDAR TACTGTAG"
motifs11="AATMAATA AGACGCAG CGGYTTGA"
motifs12="CTKCAA GCGCGTCA GTGTGTGY"
motifs13="TACCCCKA TACCTTGA"

motifs14="GAGG AGAA GATC"

TEX="arrange_plot_kinetics_correlation.extreme_ipd.celegans_v4.tex"
echo "" > $TEX
write_header >> $TEX
# High IPD motifs
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs1
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs2
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs3
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs4
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs5
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs6
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs7
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs8
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs9
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs10
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs11
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs12
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs13
echo '\end{figure}' >> $TEX
echo '\clearpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX $motifs14
echo '\end{figure}' >> $TEX
write_footer >> $TEX

lualatex $TEX
pdfcrop --hires "$(basename $TEX .tex).pdf"
