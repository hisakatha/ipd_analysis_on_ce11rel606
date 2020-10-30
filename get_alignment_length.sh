#!/usr/bin/env bash
# Input should be in SAM format without headers
cat $1 | cut -f 6 | sed -E -e "s/[0-9]+[ISHP]//g" -e "s/([0-9]+)[MDN=X]/+\1/g" -e "s/^/0/" | bc
