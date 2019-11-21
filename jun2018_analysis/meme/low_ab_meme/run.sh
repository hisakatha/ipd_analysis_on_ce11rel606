#!/usr/bin/env bash
echo "Low ab (local installation)" > description
~/memesuite/bin/meme-chip -oc . -time 300 -ccut 100 -fdesc description -order 1 -db ~/memesuite/db/motif_databases/WORM/uniprobe_worm.meme -meme-mod anr -meme-minw 4 -meme-maxw 30 -meme-nmotifs 8 -meme-searchsize 100000 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 ../../extreme_ipd.low.ab.coverage25.slop20.gff.fa
