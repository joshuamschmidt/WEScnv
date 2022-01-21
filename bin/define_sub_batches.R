#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
pacman::p_load(data.table,FNN,stringdist,cluster,svglite,ggplot2,ggrepel,factoextra,install=T,update=F)


fpkm_file=args[1]

# k-medoids
fpkm <- fread(fpkm_file,header=T)
setnames(fpkm, old=c("CHR","START","END"), new=c("chromosome","start","end"))
fpkm <- fpkm[chromosome %in% paste0("chr",1:22)]
# 20% of data for cluster inference
sample_N = round(fpkm[,.N]/5)
fpkm_sub <- fpkm[sample(.N,size=sample_N,replace=F)][,-(1:9)]
fpkm_subM <- t(as.matrix(fpkm_sub))
no_data <- which(colSums(fpkm_subM)==0)
fpkm_subM <- fpkm_subM[,-no_data]
fpkm_subM <- scale(fpkm_subM,center=T,scale=T)
k_stats <- fviz_nbclust(fpkm_subM, pam,method = 'silhouette',k.max = 10)

best_k <- k_stats$data$clusters[which(k_stats$data$y==max(k_stats$data$y))]
pam.res <- pam(fpkm_subM, best_k)
clusterDT <- data.table(id=names(pam.res$clustering),cluster=pam.res$clustering)


# k medoids plots
svglite("optimal-k-medoid-clusters.svg", width = 10, height = 10)
k_stats
dev.off()

svglite("k-medoid-clusters.svg", width = 15, height = 15)
fviz_cluster(pam.res,
             ellipse.type = "t", # Concentration ellipse
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_bw(),
             labelsize = 5,
             pointsize=0.75

)
dev.off()

# Get 100-nearest neighbors for each sample

k.param <- 100
knns <- get.knn(fpkm_subM,k=k.param,algorithm="kd_tree")

# Generate a single file for each sample listing its k-nearest neighbor sample IDs
#fpkm_subM <- cbind(fpkm_subM,rep.int(0,times=dim(fpkm_subM)[1]))
fpkm_subM <- as.data.table(fpkm_subM,keep.rownames = "id")
fpkm_cols <- grep('^V',names(fpkm_subM))

for (i in 1:fpkm_subM[,uniqueN(id)]) {
    i_sample <- fpkm_subM[i,id]
    med_cluster <- clusterDT[id==i_sample,cluster]
    nn.indexs <- knns$nn.index[i,]
    nn.clusters <- clusterDT[nn.indexs, cluster]
    # in same k-medoid cluster
    correct_cluster <- which(nn.clusters==med_cluster)
    nn.indexs <- nn.indexs[correct_cluster]
    nn.sampleids <- fpkm_subM[ nn.indexs, id ]
    # catch any very closely named - either dups or family members
    non_rel <- which(stringsim(i_sample, nn.sampleids,method = 'lcs') < 0.85)
    nn.indexs <- nn.indexs[non_rel]
    nn.sampleids <- nn.sampleids[non_rel]
    fname <- paste(i_sample, ".", k.param, "nns.txt", sep="")
    fwrite(x=data.table(sample_id=c(i_sample,nn.sampleids)),
         file=fname,col.names=T, row.names=F,quote=F,sep="\t")
    fpkm_subM[i,n_ref:=length(nn.sampleids)]
    # mean distance of reference samples
    center <- colMeans(fpkm_subM[nn.indexs, ..fpkm_cols]);
    mean_dist <- as.numeric(dist(rbind(as.numeric(fpkm_subM[i, ..fpkm_cols]), as.numeric(center))))
    fpkm_subM[i,DistanceToReferenceSet:=mean_dist]
}
e = ecdf(fpkm_subM$DistanceToReferenceSet)
fpkm_subM[,`ecdf(NN_distance)`:=e(DistanceToReferenceSet)]
clusterDT <- clusterDT[fpkm_subM[, !..fpkm_cols],on=.(id)]
clusterDT[,`k-medoid Cluster`:=paste0("cluster_",cluster)]
clusterDT[,`top 10% ref distance`:=ifelse(`ecdf(NN_distance)`> 0.9,"YES","NO")]
clusterDT[,label:=ifelse(`ecdf(NN_distance)`> 0.9,id,"")]
# Plot distance distribution

svglite("knn-distance-cluster-ref_n.svg", width = 15, height = 10)
ggplot(clusterDT, aes(x=DistanceToReferenceSet,y=`k-medoid Cluster`,color=`k-medoid Cluster`)) + geom_jitter(aes(size = n_ref),alpha=0.5) + theme_bw()
dev.off()


svglite("knn-distance-distribution.svg", width = 10, height = 10)
ggplot(clusterDT, aes(x=DistanceToReferenceSet,y=`ecdf(NN_distance)`,color=`k-medoid Cluster`)) + geom_point(alpha=0.5) + theme_bw() + geom_text_repel(data = clusterDT, aes(DistanceToReferenceSet,`ecdf(NN_distance)`, label=label ),max.overlaps = 50,size=2)
dev.off()

# write table of metrics.....
fwrite(x=clusterDT[,.(sample_id=id,cluster,n_ref,DistanceToReferenceSet,`ecdf(NN_distance)`)],
       file="knn-kmedoid-stats.txt",col.names=T, row.names=F,quote=F,sep="\t")
