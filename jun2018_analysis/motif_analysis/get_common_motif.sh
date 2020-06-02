#!/usr/bin/env bash
cat ../meme/high_ab_celegans/summary.sort.txt ../meme/high_cd_celegans/summary.sort.txt | grep -v CONSENSUS | grep -v -e "^$" | sort | uniq -d > common_motif.high_celegans.txt
cat ../meme/low_ab_celegans/summary.sort.txt ../meme/low_cd_celegans/summary.sort.txt | grep -v CONSENSUS | grep -v -e "^$" | sort | uniq -d > common_motif.low_celegans.txt
