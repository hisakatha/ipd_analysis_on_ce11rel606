#!/usr/bin/env bash
awk '$0 ~ /^#/{print} $3 == "m6A"{cov = gensub(/.*coverage=([0-9]+).*/, "\\1", 1, $9); if(cov + 0 >= 25){print}}' "$@"
