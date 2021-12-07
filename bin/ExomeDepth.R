#!/usr/bin/env Rscript 
args = commandArgs(trailingOnly=TRUE)
library(stringdist)
library(data.table)
library(ExomeDepth)

# from ARGV
test_sample <- args[1]
counts <- fread(args[2],header=T)


counts[,width:=end-start]
# construct reference set for the sample....
exomes <- names(counts)[!grepl("user|chr|start|end",names(counts))]
test_sample <- exomes[1]
ref_samples <- exomes[!exomes %in% test_sample]
test.exome <- counts[,get(test_sample)]
ref.exomes <- as.matrix(counts[, ..ref_samples])

my.choice <- select.reference.set(test.counts = test.exome,
                                  reference.counts = ref.exomes,
                                  bin.length = (counts$end - counts$start)/1000,
                                  data = counts[,.(MAP=ifelse(user_5 %in% NaN,0,user_5),GC=user_7,width)],
                                  formula = 'cbind(test, reference) ~ 1 + MAP + GC + width',
                                  phi.bins = 1)
best_ref_set <- my.choice$reference.choice
my.matrix <- as.matrix(counts[, ..best_ref_set])
my.reference.selected <- apply(X = my.matrix,
                               MAR = 1,
                               FUN = sum)

all.exons <- new('ExomeDepth',
                 test = test.exome,
                 reference = my.reference.selected,
                 data = counts[,.(MAP=ifelse(user_5 %in% NaN,0,user_5),GC=user_7)],
                 formula = 'cbind(test, reference) ~ 1 + MAP + GC',
                 phi.bins=1)

counts$CHR <- gsub(as.character(counts$chr),
                                   pattern = 'chr',
                                   replacement = '')


all.exons <- CallCNVs(x = all.exons,
                      transition.probability = 10^-4,
                      chromosome = counts$CHR,
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
all.exons@CNV.calls[3:6,]


stringdist("abc","abcd", method = "lv")
