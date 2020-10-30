#!/usr/bin/env bash
high_motifs="AAAABT AGAGAGTA AGCTATAT AGGCAGGC ATGGGAYA ATTGTTAC CAAACTAC CAGYTG CCAATCAG CRACGAS DCGAGACC DGCTTC DHNDATCAGCTGATSD DNKAMCAGSTGKYMNN GAAGGATC GATATRGY GCCVCGCC GCGACCTA GCGCGCGC GGHGGY GTAGATCA GTATCGTA HSDAACAKCTGTCHHW KVGACAGCTGTCRDVV MRWWNCACCTGCTNDNN NBMAACATMTGGTNKM NDNGRCACGTGTCNNND NDRHYACGCRTGKYHND NMDAACAGCTGTYNHN NNAACABCTGTCVNNA NNDACCACGTGGTYND RKRCGCGTGTCCCWN RNGRCGCGTGTCCYND RNMVCASCTGBGBCMT TGACGTCA TGGTGSA TNDRSCACGTGSBANN VRKRMCACCTGCTNVHH
WNAHGGCACGTGYCNDA WNDATCACGTGAYHWDN YBKRVARWAAADWG YGGAR"
low_motifs="AATMAATA ADNRATWAGWANWH AGACGCAG AWBWWTWAAAANWA CGGYTTGA CTGDAR CTKCAA CWHAAAAAAWAAAA GCGCGTCA GTGTGTGY HSDAACAKCTGTCHHW KVGACAGCTGTCRDVV NBMAACATMTGGTNKM NBNAATATAAWVHT NMDAACAGCTGTYNHN NNAACABCTGTCVNNA RNMVCASCTGBGBCMT TACCCCKA TACCTTGA TGACGTCA WAYYAWATWTRKVV WHWAATAAAACNVW YBKRVARWAAADWG"
motifs="$high_motifs $low_motifs"
for motif in $motifs
do
    if [[ -f $motif/PATTERN ]]; then
        echo "Skipping making $motif/PATTERN"
    else
        mkdir -p $motif
        echo $motif > $motif/PATTERN
    fi
done
