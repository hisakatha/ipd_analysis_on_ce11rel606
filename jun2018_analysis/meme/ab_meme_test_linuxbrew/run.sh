#!/usr/bin/env bash
MEMEPREFIX=/home/hisakatha/.linuxbrew/opt/meme/libexec
$MEMEPREFIX/libexec/meme-5.0.4/fasta-get-markov -dna -m 1 ../../extreme_ipd.low.ab.coverage25.slop20.gff.fa ./background

#alias sh=bash
#$MEMEPREFIX/bin/meme ../../extreme_ipd.low.ab.coverage25.slop20.gff.fa -p 8 -oc meme_out -mod anr -nmotifs 8 -minw 4 -maxw 30 -bfile ./background -dna -searchsize 100000 -revcomp
$MEMEPREFIX/bin/meme ../../extreme_ipd.low.ab.coverage25.slop20.gff.fa -p 20 -oc meme_out -mod anr -nmotifs 3 -bfile ./background -dna -searchsize 100000 -revcomp
