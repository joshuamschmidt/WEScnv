# GATK

Make a singularity image based on docker hub:


`singularity build /share/ClusterShare/software/contrib/jossch/singularity_images/gatk:4.2.3.0  docker://broadinstitute/gatk:4.2.3.0`

`singularity build /share/ClusterShare/software/contrib/jossch/singularity_images/picard:2.26.6 docker://broadinstitute/picard:2.26.6`


#!/bin/bash
echo "Begin ${task_file}"

# depends
module load samtools
module load bedtools
module load R/4.0.2-foss-2020a
picard_dir="/home/schm0229/picard/build/libs"

# env vars
export INSERT_SIZE=200
export CLAMMS_DIR="/home/schm0229/clamms"

# tmp workspace
tmp_dir="/local/$SLURM_JOBID/$SLURM_ARRAY_TASK_ID"
mkdir -p ${tmp_dir}

# data + results folders
task_file=${1}
b_name=$(basename -s ".bam" $task_file)

bams="/scratch/user/schm0229/glaucoma_exomes/bam"
targets="/scratch/user/schm0229/cnv_calling/targets"
masks="/scratch/user/schm0229/cnv_calling/genome_masks"
ref="/scratch/user/schm0229/rna_seq/reference"
CLAMMS_WORK="/scratch/user/schm0229/cnv_calling/CLAMMS"
mkdir -p $CLAMMS_WORK

##--- pre
#sed -e 's/^chr//g' $targets/Exome-Agilent_V4.bed | sort -u -k1,1 -k2,2 -k3,3 | sort -k1,1 -k2,2n > $targets/Exome-Agilent_V4-CLAMMS.bed
#sed -e 's/^>chr/>/g' $ref/hg19.fa  > $ref/hg19-nochr.fa

#$CLAMMS_DIR/annotate_windows.sh \
#${targets}/Exome-Agilent_V4-CLAMMS.bed \
#${ref}/hg19-nochr.fa    \
#${masks}/CLAMMSmappability.sort.bed \
#$INSERT_SIZE \
#$CLAMMS_DIR/data/clamms_special_regions.hg19.bed > $CLAMMS_WORK/Exome-Agilent_V4-CLAMMS-windows.bed
###-----

cp ${bams}/${task_file}* ${tmp_dir}
cp ${bams}/${task_file}.bai ${tmp_dir}
cd ${tmp_dir}


#---- get the required picard metrics
metrics_out=$CLAMMS_WORK/picard_qc
mkdir -p ${metrics_out}

# craziness with bam headers and picard tools. need to do this.
samtools view -H ${b_name}.bam | grep SQ | tr ':' '\t' | cut -f3 > ${tmp_dir}/order.txt
#samtools faidx ${ref}/hg19.fa $(cat ${tmp_dir}/order.txt) > ${tmp_dir}/hg19.bam-order.fa
samtools faidx ${ref}/broad.25chr.genome.fa $(cat ${tmp_dir}/order.txt) > ${tmp_dir}/hg19.bam-order.fa
samtools faidx ${tmp_dir}/hg19.bam-order.fa

java -Xmx12g -jar ${picard_dir}/picard.jar CreateSequenceDictionary \
R=${tmp_dir}/hg19.bam-order.fa \
O=${tmp_dir}/hg19.bam-order.fa.dict

java -Xmx12g -jar ${picard_dir}/picard.jar BedToIntervalList \
I=${targets}/Exome-Agilent_V4.bed \
O=${tmp_dir}/Exome-Agilent_V4.interval_list \
SD=${tmp_dir}/hg19.bam-order.fa.dict

#---- collect HSmetrics
java -Xmx12g -jar ${picard_dir}/picard.jar CollectHsMetrics \
I=${b_name}.bam \
O=${b_name}_HSmetrics.txt \
R=${tmp_dir}/hg19.bam-order.fa \
BAIT_INTERVALS=${tmp_dir}/Exome-Agilent_V4.interval_list  \
TARGET_INTERVALS=${tmp_dir}/Exome-Agilent_V4.interval_list

#---- collect InsertSizeMetrics
java -Xmx12g -jar ${picard_dir}/picard.jar CollectInsertSizeMetrics \
I=${b_name}.bam \
O=${b_name}_InsertSizemetrics.txt \
H=${b_name}_InsertSizeHist.pdf


#---- sample coverage + normalisation
samtools bedcov -Q 30 \
<(awk '{print "chr"$0}' $CLAMMS_WORK/Exome-Agilent_V4-CLAMMS-windows.bed) \
$tmp_dir/${b_name}.bam \
| awk '{ printf "%s\t%d\t%d\t%.6g\n", $1, $2, $3, $NF/($3-$2); }' \
| sed -e 's/^chr//g' > $tmp_dir/${b_name}.coverage.bed

$CLAMMS_DIR/normalize_coverage \
$tmp_dir/${b_name}.coverage.bed \
$CLAMMS_WORK/Exome-Agilent_V4-CLAMMS-windows.bed \
> $tmp_dir/${b_name}.norm.cov.bed

#---- transfer results to scratch
cp $tmp_dir/${b_name}*.bed $CLAMMS_WORK/CLAMMs_cov
cp $tmp_dir/${b_name}_InsertSize* ${metrics_out}
cp $tmp_dir/${b_name}_HSmetrics.txt ${metrics_out}

#--- tidy up
rm -rf ${tmp_dir}




##

#!/usr/bin/env bash

# need to set env variables
# export SCRIPTS=//
# export METRICS_DIR=//
# export nnK=10
files=("$@")
if ((${#files[@]} == 0)); then
  echo "error: no files to work on!";
  exit 1
fi

TODAY=$(date +"%Y-%m-%d")

echo -e "SAMPLE\tONBAITVSSELECTED\tPCTPFUQREADS\tPCTTARGETBASES1X\tPCTTARGETBASES10X\tPCTTARGETBASES50X\tATDROPOUT\tGCDROPOUT\tMEANINSERTSIZE" > $METRICS_DIR/${TODAY}_picardQC.txt

for file in "${files[@]}"; do
  sample=$(basename -s "_HSmetrics.txt" $file);
  #echo $sample;
  HS=($(sed '8q;d' ${sample}_HSmetrics.txt | cut -f9,32,46,48,52,63,64));
  IS=$(sed '8q;d' ${sample}_InsertSizemetrics.txt | cut -f6 );
  OUT=($sample "${HS[@]}" $IS);
  echo "${OUT[*]}" | tr ' ' '\t' >> $METRICS_DIR/${TODAY}_picardQC.txt;
done

Rscript $SCRIPTS/cluster_metrics.Rscript $METRICS_DIR/${TODAY}_picardQC.txt $nnK


##

#!/usr/bin/env RScript
# ARGV
args <- commandArgs(trailingOnly = TRUE)
input=args[1]
k.param=as.integer(args[2])

# DEPENDS
library('FNN')
library('ggplot2')
library('data.table')
setDTthreads(threads=2)
library('ggrepel')
library('svglite')

# WORK
qcs <- fread(input,header=T)
out_prefix <- strsplit(x = input,split = '-p')[[1]][1]
# Create a scaled copy of the data frame
qcs.scaled <- qcs
cols <- names(qcs.scaled)[2:ncol(qcs.scaled)]
# any columns with no variation in metric?
cond <- qcs.scaled[, lapply(.SD, function(x) uniqueN(x)==1), .SDcols = cols]
if(any(cond==TRUE)){
  qcs.scaled[, which(cond==TRUE)+1 := NULL][]
  cols <- names(qcs.scaled)[2:ncol(qcs.scaled)]
}
# scale data
for (j in cols) set(qcs.scaled, j = j, value = scale(qcs.scaled[[j]]))

# Get k-nearest neighbors for each sample
knns <- get.knn(qcs.scaled[,2:ncol(qcs.scaled)],k=k.param,algorithm="kd_tree")

# Generate a single file for each sample listing its k-nearest neighbor sample IDs
for (i in 1:nrow(qcs.scaled)) {
    fname <- paste(qcs.scaled$SAMPLE[i], ".", k.param, "nns.txt", sep="")
    nn.sampleids <- qcs.scaled$SAMPLE[ knns$nn.index[i,] ]
    write.table(nn.sampleids, fname, quote=F, row.names=F, col.names=F)
}


# To check how well each sample's kNNs fit, compute the distance to its kNN cluster mean
qcs.scaled$DistanceToClusterMean <- sapply(1:nrow(qcs.scaled),  function(x) {
    this.knns <- knns$nn.index[x,];
    center <- colMeans(qcs.scaled[this.knns, 2:ncol(qcs.scaled)]);
    return(as.numeric(dist(rbind(as.numeric(qcs.scaled[x, 2:ncol(qcs.scaled)]), as.numeric(center)))))
})
qcs.scaled$label <- qcs.scaled$SAMPLE
qcs.scaled[DistanceToClusterMean < 1,label:=""]

e = ecdf(qcs.scaled$DistanceToClusterMean)
qcs.scaled[,ed:=e(DistanceToClusterMean)]
# plot the result
p <- ggplot(qcs.scaled, aes(DistanceToClusterMean)) +
  stat_ecdf() + theme_classic() + geom_text_repel(data = qcs.scaled, aes(DistanceToClusterMean,ed, label=label ),max.overlaps = 50,size=2)





# Plot distance distribution
ggsave(filename=paste0(out_prefix,"_nnK-",k.param,"_qc.svg"),
    plot= p,
    device="svg")


# save table of metrics with distance to cluster.
fwrite(x=qcs.scaled,
    file=paste0(out_prefix,"_nnK-",k.param,"_picardQC.scaled.txt"),
    col.names=T,
    row.names=F,
    quote=F,
    sep='\t')

