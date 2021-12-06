#!/usr/bin/env bash
infile=${1}
prefix=${2}

# params for hmm. change Pr enter CNV to 1e-06 (1e-08)
echo -e "1e-06\t6\t70\t-3\t1\t0\t1\t3\t1" > params.txt

n_cols=$(zcat < "$infile" |head -n1 | cut -f 10- | awk '{print NF}')
n_rows=$(zcat < "$infile" | sed '1d'  | wc -l)
zcat < "$infile" | awk -v chr=$chr 'BEGIN {OFS='\t'; print "GATK._mean_cvg";} {if ($5 >= 0.2 && $5 <= 0.8) print $1":"$2"-"$3}' | sed '/^$/d' | transpose -t --limit "$n_rows"x1 | sed -e 's/[[:space:]]*$//g' > tmp_header.txt
zcat < "$infile" | awk -v chr=$chr 'NR==1 || ($5 >= 0.2 && $5 <= 0.8)' | cut -f 10-  | transpose -t --limit "$n_rows"x"$n_cols" | sed -e 's/[[:space:]]*$//g' > tmp_data.txt

cat tmp_header.txt tmp_data.txt > "$prefix".txt
rm tmp_header.txt && rm tmp_data.txt

#docker run -u $(id -u):$(id -g) --rm -v $(pwd):/app --rm -ti joshmschmidt/xhmm /bin/bash

DATA="$prefix"

xhmm --matrix -r ./"$DATA".txt --centerData --centerType target \
-o ./"$DATA".filtered_centered.txt \
--outputExcludedTargets ./"$DATA".filtered_centered.txt.filtered_targets.txt \
--outputExcludedSamples ./"$DATA".filtered_centered.txt.filtered_samples.txt \
--minTargetSize 10 --maxTargetSize 10000 \
--minMeanTargetRD 10 --maxMeanTargetRD 500 \
--minMeanSampleRD 20 --maxMeanSampleRD 200 \
--maxSdSampleRD 150

xhmm --PCA -r ./"$DATA".filtered_centered.txt --PCAfiles ./"$DATA".RD_PCA

xhmm --normalize -r ./"$DATA".filtered_centered.txt --PCAfiles ./"$DATA".RD_PCA \
--normalizeOutput ./"$DATA".PCA_normalized.txt \
--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7

xhmm --matrix -r ./"$DATA".PCA_normalized.txt --centerData --centerType sample --zScoreData \
-o ./"$DATA".PCA_normalized.filtered.sample_zscores.txt \
--outputExcludedTargets ./"$DATA".PCA_normalized.filtered.sample_zscores.txt.filtered_targets.txt \
--outputExcludedSamples ./"$DATA".PCA_normalized.filtered.sample_zscores.txt.filtered_samples.txt \
--maxSdTargetRD 30

xhmm --matrix -r ./"$DATA".txt \
--excludeTargets ./"$DATA".filtered_centered.txt.filtered_targets.txt \
--excludeTargets ./"$DATA".PCA_normalized.filtered.sample_zscores.txt.filtered_targets.txt \
--excludeSamples ./"$DATA".filtered_centered.txt.filtered_samples.txt \
--excludeSamples ./"$DATA".PCA_normalized.filtered.sample_zscores.txt.filtered_samples.txt \
-o ./"$DATA".same_filtered.txt

xhmm --discover -p ./params.txt \
-r ./"$DATA".PCA_normalized.filtered.sample_zscores.txt -R ./"$DATA".same_filtered.txt \
-c ./"$DATA".xcnv -a ./"$DATA".aux_xcnv -s ./"$DATA"

xhmm --genotype -p ./params.txt \
-r ./"$DATA".PCA_normalized.filtered.sample_zscores.txt -R ./"$DATA".same_filtered.txt \
-g ./"$DATA".xcnv \
-v ./"$DATA".xhmm.calls.vcf
