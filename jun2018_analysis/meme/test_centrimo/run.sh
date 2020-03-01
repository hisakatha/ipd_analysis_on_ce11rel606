#!/usr/bin/env bash
$HOME/memesuite/bin/centrimo -seqlen 0 -verbosity 1 -oc centrimo_out -bfile ../low_ab_meme/background -score 5.0 -ethresh 10.0 ../low_ab_meme/extreme_ipd.low.ab.coverage25.slop20.gff.fa dreme.xml
