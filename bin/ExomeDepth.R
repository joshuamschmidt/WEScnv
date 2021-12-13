#!/usr/bin/env Rscript 
args = commandArgs(trailingOnly=TRUE)
library(stringdist)
library(data.table)
library(ExomeDepth)

# from ARGV
counts <- fread(args[1],header=T)
test_sample <- args[2]
candidate_ref_set <- args[3]

counts <- counts[chr %in% paste0("chr",1:22)]
# construct reference set for the sample....
exomes <- names(counts)[!grepl("user|chr|CHR|start|end",names(counts))]
exome_matrix <- as.matrix(counts[, ..exomes])
colnames(exome_matrix) <- exomes
rownames(exome_matrix) <- counts$user_4
exome_matrix <- t(exome_matrix)
exome_matrix <- exome_matrix / rowSums(exome_matrix)

ref_samples <- exomes[!exomes %in% test_sample]
test.exome <- counts[,get(test_sample)]
ref.exomes <- as.matrix(counts[, ..ref_samples])
my.choice <- select.reference.set(test.counts = test.exome,
                                  reference.counts = ref.exomes,
                                  bin.length = (counts$end - counts$start)/1000,
                                  data = counts[,.(MAP=ifelse(user_5 %in% NaN,0,user_5),GC=user_6,GC500=user_7,width=user_8,merged=user_9)],
                                  formula = 'cbind(test, reference) ~ 1 + MAP + GC + GC500 + width + merged',
                                  phi.bins = 1,
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
                 data = counts[,.(MAP=ifelse(user_5 %in% NaN,0,user_5),GC=user_6,GC500=user_7,width=user_8,merged=user_9)],
                 formula = 'cbind(test, reference) ~ 1 + MAP + GC + GC500 + width + merged',
                 phi.bins=1)


# counts$CHR <- gsub(as.character(counts$chr),
#                                    pattern = 'chr',
#                                   replacement = '')


all.exons <- CallCNVs(x = all.exons,
                      transition.probability = 10^-4,
                      chromosome = gsub(as.character(counts$chr),
                                        pattern = 'chr',
                                        replacement = ''),
                      start = counts$start,
                      end = counts$end,
                      name = counts$user_4)



exons.hg38.GRanges <- GenomicRanges::GRanges(seqnames = counts$CHR,
                                             IRanges::IRanges(start=counts$start,end=counts$end),
                                             names = counts$user_4)
all.exons <- AnnotateExtra(x = all.exons,
                           reference.annotation = exons.hg38.GRanges,
                           min.overlap = 0.0001,
                           column.name = 'exons.hg38')
calls <- data.table(all.exons@CNV.calls[ order ( all.exons@CNV.calls$BF, decreasing = TRUE),])[order(-BF)]


stringdist("abc","abcd", method = "lv")
