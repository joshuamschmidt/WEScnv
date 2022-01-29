#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
pacman::p_load('data.table','ExomeDepth',install = F, update = F)
# from ARGV
counts <- fread(args[1],header=T)
setnames(counts, old=c("CHR","START","END"), new=c("chromosome","start","end"))
counts <- counts[chromosome == "X"]
counts <- counts[!MAP100 %in% NaN]
counts <- counts[order(chromosome,start,end)]
# test and reference exomes
test_candidate_ref_set <- fread(args[2],header=T)
test_sample <- test_candidate_ref_set[1,sample_id]

# construct reference set for the sample....
candidate_ref_set <- test_candidate_ref_set[!sample_id %in% test_sample,sample_id]
# bioSex of test and references
bio_sex <- fread(args[3], header=T)
test_sex <- bio_sex[Sample_ID==test_sample,bioSex]
candidate_ref_set <- bio_sex[bioSex==test_sex][Sample_ID %in% candidate_ref_set,Sample_ID]

test.exome <- counts[,get(test_sample)]
ref.exomes <- as.matrix(counts[, ..candidate_ref_set])
my.choice <- select.reference.set(test.counts = test.exome,
                                  reference.counts = ref.exomes,
                                  bin.length = (counts$end - counts$start)/1000,
                                  data = counts[,.(MAP100,GC,GC500,WIDTH)],
                                  formula = 'cbind(test, reference) ~ 1 + MAP100 + GC+ GC500 + WIDTH',
                                  n.bins.reduced = 5000)

best_ref_set <- my.choice$reference.choice
n_ref_set <- length(best_ref_set)
my.matrix <- as.matrix(counts[, ..best_ref_set])
my.reference.selected <- apply(X = my.matrix,
                               MAR = 1,
                               FUN = sum)

all.exons <- new('ExomeDepth',
                 test = test.exome,
                 reference = my.reference.selected,
                 data = counts[,.(MAP100,GC,GC500,WIDTH)],
                 formula = 'cbind(test, reference) ~ 1 + MAP100 + GC + GC500 + WIDTH',
                 phi.bins=1)


all.exons <- CallCNVs(x = all.exons,
                      transition.probability = 5e-04,
                      chromosome = gsub(as.character(counts$chr),
                                        pattern = 'chr',
                                        replacement = ''),
                      start = counts$start,
                      end = counts$end,
                      name = counts$TARGET)



exons.hg38.GRanges <- GenomicRanges::GRanges(seqnames = gsub(as.character(counts$chromosome),
                                                             pattern = 'chr',
                                                             replacement = '') ,
                                             IRanges::IRanges(start=counts$start,end=counts$end),
                                             names = counts$TARGET)
all.exons <- AnnotateExtra(x = all.exons,
                           reference.annotation = exons.hg38.GRanges,
                           min.overlap = 0.00001,
                           column.name = 'exons.hg38')
calls <- data.table(all.exons@CNV.calls[ order ( all.exons@CNV.calls$BF, decreasing = TRUE),])[order(-BF)]
calls <- calls[BF >0]
reorder <- names(calls)[c(7,5:6,3,8:13,1,2,4)]
setcolorder(calls, reorder)
calls[,start:=start-1]
calls[,chromosome:=paste0("chr",chromosome)]
setorder(calls, chromosome, start, end)
outfile <- paste(test_sample,"ExomeDepth-CNV.calls.bed",sep="_")
fwrite(calls,outfile,sep="\t",col.names = F, row.names = F, quote = F)
outfile <- paste(test_sample,"ExomeDepth-CNV.reference.txt",sep="_")
write.table(x=c(test_sex,best_ref_set),file=outfile,sep="\n",col.names = F, row.names = F, quote = F)
