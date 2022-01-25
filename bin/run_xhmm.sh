#!/usr/bin/env ash
coverageFile=${1}

cluster=$(echo $prefix | tr '_' '\t' | cut -f2)
prefix=$(echo 'cluster_'"$cluster"'_xhmm')
# params for hmm. change Pr enter CNV to 1e-07 (1e-08)
echo -e "1e-07\t6\t70\t-3\t1\t0\t1\t3\t1" > params.txt

#docker run -u $(id -u):$(id -g) --rm -v $(pwd):/app --rm -ti joshmschmidt/xhmm /bin/bash

xhmm --matrix -r ./"$coverageFile" --centerprefix --centerType target \
-o ./"$prefix".filtered_centered.txt \
--outputExcludedTargets ./"$prefix".filtered_centered.txt.filtered_targets.txt \
--outputExcludedSamples ./"$prefix".filtered_centered.txt.filtered_samples.txt \
--minTargetSize 10 --maxTargetSize 10000 \
--minMeanTargetRD 15 --maxMeanTargetRD 1000 \
--minMeanSampleRD 15 --maxMeanSampleRD 200 \
--maxSdSampleRD 150

xhmm --PCA -r ./"$prefix".filtered_centered.txt --PCAfiles ./"$prefix".PCA

xhmm --normalize -r ./"$prefix".filtered_centered.txt --PCAfiles ./"$prefix".PCA \
--normalizeOutput ./"$prefix".PCA_normalized.txt \
--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7

xhmm --matrix -r ./"$prefix".PCA_normalized.txt --centerprefix --centerType sample --zScoreprefix \
-o ./"$prefix".PCA_normalized.filtered.sample_zscores.txt \
--outputExcludedTargets ./"$prefix".PCA_normalized.filtered.sample_zscores.txt.filtered_targets.txt \
--outputExcludedSamples ./"$prefix".PCA_normalized.filtered.sample_zscores.txt.filtered_samples.txt \
--maxSdTargetRD 30

xhmm --matrix -r ./"$prefix".txt \
--excludeTargets ./"$prefix".filtered_centered.txt.filtered_targets.txt \
--excludeTargets ./"$prefix".PCA_normalized.filtered.sample_zscores.txt.filtered_targets.txt \
--excludeSamples ./"$prefix".filtered_centered.txt.filtered_samples.txt \
--excludeSamples ./"$prefix".PCA_normalized.filtered.sample_zscores.txt.filtered_samples.txt \
-o ./"$prefix".same_filtered.txt

xhmm --discover -p ./params.txt \
-r ./"$prefix".PCA_normalized.filtered.sample_zscores.txt -R ./"$prefix".same_filtered.txt \
-c ./"$prefix".xcnv -a ./"$prefix".aux_xcnv -s ./"$prefix"

xhmm --genotype -p ./params.txt \
-r ./"$prefix".PCA_normalized.filtered.sample_zscores.txt -R ./"$prefix".same_filtered.txt \
-g ./"$prefix".xcnv \
-v ./"$prefix".xhmm.calls.vcf

