#!/usr/bin/env bash
high_motifs="AAAABT AGAGAGTA AGCTATAT AGGCAGGC ATGGGAYA ATTGTTAC CAAACTAC CAGYTG CCAATCAG CRACGAS DCGAGACC DGCTTC DHNDATCAGCTGATSD DNKAMCAGSTGKYMNN GAAGGATC GATATRGY GCCVCGCC GCGACCTA GCGCGCGC GGHGGY GTAGATCA GTATCGTA HSDAACAKCTGTCHHW KVGACAGCTGTCRDVV MRWWNCACCTGCTNDNN NBMAACATMTGGTNKM NDNGRCACGTGTCNNND NDRHYACGCRTGKYHND NMDAACAGCTGTYNHN NNAACABCTGTCVNNA NNDACCACGTGGTYND RKRCGCGTGTCCCWN RNGRCGCGTGTCCYND RNMVCASCTGBGBCMT TGACGTCA TGGTGSA TNDRSCACGTGSBANN VRKRMCACCTGCTNVHH
WNAHGGCACGTGYCNDA WNDATCACGTGAYHWDN YBKRVARWAAADWG YGGAR"
low_motifs="AATMAATA ADNRATWAGWANWH AGACGCAG AWBWWTWAAAANWA CGGYTTGA CTGDAR CTKCAA CWHAAAAAAWAAAA GCGCGTCA GTGTGTGY HSDAACAKCTGTCHHW KVGACAGCTGTCRDVV NBMAACATMTGGTNKM NBNAATATAAWVHT NMDAACAGCTGTYNHN NNAACABCTGTCVNNA RNMVCASCTGBGBCMT TACCCCKA TACCTTGA TGACGTCA WAYYAWATWTRKVV WHWAATAAAACNVW YBKRVARWAAADWG"

motifs_high2="ACGCRTG ATCAGCTG RGTA"
motifs_low2="ACATMTGG CTGDAR TACTGTAG"
motifs_select1="GAGG BGAGG AGAA ASAAD"
motifs_select2="GATC"
motifs_select3="ANNNATCAGCTG CNNNATCAGCTG GNNNATCAGCTG TNNNATCAGCTG ATCAGCTGATCAGCTG"
motifs_select4="ANNGATC CNNGATC GNNGATC TNNGATC GATCGATC"
# REBASE E. coli REL606
motifs_select5="ATGCAT TGANNNNNNNNTGCT"
# From Blow et al. (2016) https://doi.org/10.1371/journal.pgen.1005854
motifs_select6="ATTAAT AATT CATG RAATTY ATCGAT CTCGAG TTAA ATGCAT GANTC AGCT CTAG"

old_targets="ATWTTB CGRAACCCG CTGRAA GAGD GCGC GGTCTCGC GKCGGY GTAS HKGAGCA HRGTA RGTAB RTAB TATAA TTTT YGGWR"

motifs2="$motifs_high2 $motifs_low2 $motifs_select1 $motifs_select2 $motifs_select3 $motifs_select4 $motifs_select5 $motifs_select6"

motifs="$high_motifs $low_motifs $motifs2"
for motif in $motifs
do
    if [[ -f $motif/PATTERN ]]; then
        echo "Skipping making $motif/PATTERN"
    else
        mkdir -p $motif
        echo $motif > $motif/PATTERN
    fi
done

# GGN_4
DIR=GGN_4; mkdir -p $DIR && echo GGNGGNGGNGGN > $DIR/PATTERN
DIR=GGN_10; mkdir -p $DIR && echo GGNGGNGGNGGNGGNGGNGGNGGNGGNGGN > $DIR/PATTERN
DIR=AGGCTT_4; mkdir -p $DIR && echo AGGCTTAGGCTTAGGCTTAGGCTT > $DIR/PATTERN
DIR=AGGY_7; mkdir -p $DIR && echo AGGYAGGYAGGYAGGYAGGYAGGYAGGY > $DIR/PATTERN