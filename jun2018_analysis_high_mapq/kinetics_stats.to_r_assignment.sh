#!/usr/bin/env bash
# input: kinetics_stats.csv or kinetics_stats.c_elegans.csv
awk -v FS=, 'NR>1{prefix = $1 "_" $2 "_" $3 "_"; print prefix "mean <- " $4; print prefix "var <- " $5 * $5; print prefix "size <- " $6}'
