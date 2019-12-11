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

write_include_for_motif () {
    original_pdf="$1/plot_kinetics_correlation.pdf"
    cropped_pdf="cropped_kinetics_correlation.$1.pdf"
    pdfcrop --hires "$original_pdf" "$cropped_pdf" >&2
    echo '\begin{minipage}[t]{0.33\linewidth}'
    echo '\verb'"|$1|"'\\'
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

TEX="arrange_plot_kinetics_correlation.extreme_ipd.tex"
echo "" > $TEX
write_header >> $TEX
# High IPD motifs
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX "ACGCRTG" "ATCAGCTG" "GGN_10" "GGN_4" "HKGAGCA" "HRGTA" "RGTAB"
echo '\end{figure}' >> $TEX
# Low IPD motifs
echo '\newpage' >> $TEX
echo '\begin{figure}[t]' >> $TEX
write_all_for_motifs $TEX "ACATMTGG" "AGGCTT_4" "AGGY_7" "CTGDAR" "RTAB" "TACTGTAG"
echo '\end{figure}' >> $TEX
write_footer >> $TEX

lualatex $TEX
pdfcrop --hires "$(basename $TEX .tex).pdf"
