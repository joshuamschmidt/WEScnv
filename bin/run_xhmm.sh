#!/usr/bin/env ash
countsFile=${1}
prefix=${2}
group_info=${3}
group_id=${4}

# params for hmm. change Pr enter CNV to 1e-07 (1e-08)
echo -e "1e-07\t6\t70\t-3\t1\t0\t1\t3\t1" > params.txt

# get columns we need; chr, pos, GC, MAP + samples in batch.
gzip -dc "$countsFile" | cut -f1-3,6,5 > tmp_info.tsv

keep_cols=$(awk -v group="$group_id" '$2==group {print $1}' "$group_info" | head -n 3 | tr '\n' ',')

awk -F'\t' -v cols="$keep_cols" '(NR==1){n=split(cols,cs,",");for(c=1;c<=n;c++){for(i=1;i<=NF;i++)if($(i)==cs[c])ci[c]=i}}{for(i=1;i<=n;i++)printf "%s" FS,$(ci[i]);printf "\n"}' <(gzip -dc "$countsFile") > tmp_counts.tsv




gunzip -dc "$countsFile" | awk '$1!="chrX" && $1!="chrY"' | gzip > tmp.gz

# count array dimensions
n_cols=$(zcat < tmp.gz  |head -n1 | cut -f 10- | awk '{print NF}')
n_rows=$(zcat < tmp.gz  | sed '1d'  | wc -l)

# filter for GC content > 0.2 and < 0.8, and MAPPABILITY ==1.
zcat < tmp.gz | awk -v chr=$chr 'BEGIN {OFS='\t'; print "GATK._mean_cvg";} {if ($5 >= 0.1 && $5 <= 0.9) print $1":"$2"-"$3}' | sed '/^$/d' | transpose -t --limit "$n_rows"x1 | sed -e 's/[[:space:]]*$//g' > tmp_header.txt
zcat < tmp.gz | awk -v chr=$chr 'NR==1 || ($5 >= 0.1 && $5 <= 0.9)' | cut -f 10-  | transpose -t --limit "$n_rows"x"$n_cols" | sed -e 's/[[:space:]]*$//g' > tmp_data.txt

cat tmp_header.txt tmp_data.txt > "$prefix".txt
rm tmp_header.txt && rm tmp_data.txt

#docker run -u $(id -u):$(id -g) --rm -v $(pwd):/app --rm -ti joshmschmidt/xhmm /bin/bash

DATA="$prefix"

xhmm --matrix -r ./"$DATA".txt --centerData --centerType target \
-o ./"$DATA".filtered_centered.txt \
--outputExcludedTargets ./"$DATA".filtered_centered.txt.filtered_targets.txt \
--outputExcludedSamples ./"$DATA".filtered_centered.txt.filtered_samples.txt \
--minTargetSize 10 --maxTargetSize 10000 \
--minMeanTargetRD 15 --maxMeanTargetRD 500 \
--minMeanSampleRD 20 --maxMeanSampleRD 200 \
--maxSdSampleRD 150

xhmm --PCA -r ./"$DATA".filtered_centered.txt --PCAfiles ./"$DATA".PCA

xhmm --normalize -r ./"$DATA".filtered_centered.txt --PCAfiles ./"$DATA".PCA \
--normalizeOutput ./"$DATA".PCA_normalized.txt \
--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.9

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

rm tmp.gz
