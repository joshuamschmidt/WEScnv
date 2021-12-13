#!/usr/bin/env bash
sample_id=${1}
hsMetrics=${2}
isMetrics=${3}


HS=($(sed '8q;d' "$hsMetrics" | cut -f9,32,46,48,52,63,64));
IS=$(sed '8q;d' "$isMetrics" | cut -f6);
OUT=($sample "${HS[@]}" "$IS");
echo "${OUT[*]}" | tr ' ' '\t'
