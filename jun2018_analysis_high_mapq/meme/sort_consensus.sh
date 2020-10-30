#!/usr/bin/env bash
cat result/summary.tsv | grep -v "^#" | grep -v "^$" | grep -v CONSENSUS | cut -f 5 | sort > summary.sort.txt
