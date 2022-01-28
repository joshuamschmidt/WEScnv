#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(data.table)

coverage_file=args[1]
cluster_file=args[2]

coverage <- fread(coverage_file,header=TRUE)
clusters <- fread(cluster_file,header=TRUE)

# filter data
setnames(coverage, old=c("CHR","START","END"), new=c("chromosome","start","end"))
coverage <- coverage[chromosome %in% paste0("chr",1:22)]
coverage[,id:=paste0(chromosome,":",start,"-",end)]
all_ids <- coverage$id
coverage <- coverage[!MAP100 %in% NaN]
coverage <- coverage[MAP100 == 1]
coverage <- coverage[GC >= 0.2 & GC <= 0.8]
coverage <- coverage[WIDTH >= 10/1000 & WIDTH <= 10000/1000]
# expect only ~8000 or so exons to be excluded
filtered_ids <- all_ids[!all_ids %in% coverage$id]
unique_clusters <- sort(clusters[,unique(cluster)])

for (c in unique_clusters) {
  samples <- clusters[cluster==c,sample_id]
  if(length(samples) > 250) {
    n_chunks <- ceiling(length(samples) / 250)
    chunks <- split(samples, rep_len(1:n_chunks, length(samples)))
    for (i in 1:n_chunks){
      sub_batch <- LETTERS[i]
      tmp <- as.data.table(t(as.matrix(coverage[, ..chunks[[i]]])))
      names(tmp) <- coverage$id
      tmp$mean_cvg <- chunks[[i]]
      setcolorder(tmp, c("mean_cvg", chunks[[i]]))
      outfile <- paste0("cluster_",c,sub_batch,"_xhmm.in.txt")
      fwrite(tmp,file=outfile,sep="\t",col.names=T,row.names=F,quote=F)
      outfile <- paste0("batch_",c,sub_batch,"_XHMM.samples.txt",sep="_")
      write.table(x=chunks[[i]],file=outfile,sep="\n",col.names = F, row.names = F, quote = F)
    }
  } else {
    tmp <- as.data.table(t(as.matrix(coverage[, ..samples])))
    names(tmp) <- coverage$id
    tmp$mean_cvg <- samples
    setcolorder(tmp, c("mean_cvg", coverage$id))
    outfile <- paste0("cluster_",c,"_xhmm.in.txt")
    fwrite(tmp,file=outfile,sep="\t",col.names=T,row.names=F,quote=F)
    outfile <- paste0("batch_",c,"_XHMM.samples.txt",sep="_")
    write.table(x=samples,file=outfile,sep="\n",col.names = F, row.names = F, quote = F)
  }

}
write.table(filtered_ids,file="xhmm-filteredout.targets.txt",sep="\n",col.names=F,row.names=F,quote=F)
