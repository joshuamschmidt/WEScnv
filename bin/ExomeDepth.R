#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(stringdist)
library(data.table)
library(ExomeDepth)

# from ARGV
counts <- fread(args[1],header=T)
setnames(counts, old=c("CHR","START","END"), new=c("chromosome","start","end"))
test_sample <- args[2]
candidate_ref_set <- args[3]
candidate_ref_set <- strsplit(candidate_ref_set,split = ",")[[1]]
counts <- counts[chromosome %in% paste0("chr",1:22)]
counts <- counts[!MAP100 %in% NaN]
# construct reference set for the sample....
exomes <- names(counts)[-(1:9)]
exome_matrix <- as.matrix(counts[, ..exomes])
colnames(exome_matrix) <- exomes
rownames(exome_matrix) <- counts$TARGET
exome_matrix <- t(exome_matrix)
exome_matrix <- exome_matrix / rowSums(exome_matrix)

ref_samples <- exomes[!exomes %in% test_sample]
test.exome <- counts[,get(test_sample)]
ref.exomes <- as.matrix(counts[, ..candidate_ref_set])
#all_zeros <- which(rowSums(ref.exomes)==0)
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


stringdist("abc","abcd", method = "lv")
