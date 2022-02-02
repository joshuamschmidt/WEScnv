#!/usr/bin/env bash
coverage_file=${1}
sample_id=${2}
target_type=${3}

if [ "$target_type" == "target" ]
  then outfile="$sample_id".targetcoverage.cnn
elif [ "$target_type" == "anti_target" ]
  then outfile="$sample_id".antitargetcoverage.cnn
fi

echo -e "chromosome\tstart\tend\tgene\tdepth\tlog2" > "$outfile"
gzip -dc "$coverage_file" > tmp.bed

awk 'BEGIN {OFS="\t";} {if($5==0) logdepth="-20"; else logdepth=log($5)/log(2); print $0, logdepth}' tmp.bed >> "$outfile"
