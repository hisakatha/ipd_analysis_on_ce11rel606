#!/usr/bin/env bash
motifs="AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT"
for motif in $motifs
do
    if [[ -f $motif/PATTERN ]]; then
        echo "Skipping making $motif/PATTERN"
    else
        mkdir -p $motif
        echo $motif > $motif/PATTERN
    fi
done
