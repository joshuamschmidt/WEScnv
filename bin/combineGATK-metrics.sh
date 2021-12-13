#!/usr/bin/env ash
#sample_id=${1}
hsMetrics=${1}
isMetrics=${2}


HS=($(sed '8q;d' "$hsMetrics" | cut -f9,32,46,48,52,63,64));
IS=$(sed '8q;d' "$isMetrics" | cut -f1,3,6,7);
OUT=($sample "${HS[@]}" "$IS");
echo "${OUT[*]}" | tr ' ' '\t'
