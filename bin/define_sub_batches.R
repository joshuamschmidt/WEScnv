#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(data.table)
library(FNN)
library(stringdist)
library(cluster)
library(svglite)
library(ggplot2)
library(ggrepel)
library(factoextra)
setDTthreads(threads=1)

files=args

qc <- data.table()
for (f in files){
 tmp <- fread(f)
qc <- rbindlist(list(qc,tmp))
}

# Create a scaled copy of the data frame
qc.scaled <- as.data.frame(qc)
for (i in 2:ncol(qc.scaled)) {
    mini <- min(qc.scaled[,i])
    maxi <- max(qc.scaled[,i])
    qc.scaled[,i] <- apply(qc.scaled, 1, function(row) {
    row[[i]] <- (as.numeric(row[[i]]) - mini) / (maxi - mini)
  } )
}


# k medoids
x <- as.matrix(qc[,2:12])
row.names(x) <- qc$V1
k_stats <- fviz_nbclust(x, pam,method = 'silhouette',stand=T,k.max = 15)
svglite("optimal-k.svg", width = 10, height = 10)
k_stats
dev.off()


best_k <- which(k_stats$data$y==max(k_stats$data$y))
pam.res <- pam(x, best_k)
dd <- cbind(x, cluster = pam.res$cluster)

svglite("k-medoid-clusters.svg", width = 10, height = 10)
fviz_cluster(pam.res,
             ellipse.type = "t", # Concentration ellipse
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_bw(),
             labelsize = 4,
             pointsize = 0.5
)
dev.off()


# Get k-nearest neighbors for each sample
k.param <- 100
knns <- get.knn(qc.scaled[,c(seq(2,ncol(qc.scaled)))],k=k.param,algorithm="kd_tree")


# Generate a single file for each sample listing its k-nearest neighbor sample IDs
qc.scaled$n_ref <- 0
for (i in 1:nrow(qc.scaled)) {
    i_sample <- qc.scaled$V1[i]
    nn.sampleids <- qc.scaled$V1[ knns$nn.index[i,] ]
    closest <- which( knns$nn.dist[i,] <= 0.5)
    if(length(closest) < 20) {
      nn.sampleids <- nn.sampleids[1:20]
    } else{
      nn.sampleids <- nn.sampleids[closest]
    }
    # catch any very closely named - either dups or family members
    non_rel <- which(stringsim(i_sample, nn.sampleids,method = 'lcs') < 0.85)
    fname <- paste(qc.scaled$V1[i], ".", k.param, "nns.txt", sep="")
    fwrite(x=data.table(sample_id=c(i_sample,nn.sampleids[non_rel])),
           file=fname,col.names=T, row.names=F,quote=F,sep="\t")
    qc.scaled[i,]$n_ref <- length(nn.sampleids[non_rel])
}

# plot of the cum distributon of mean differences....
qc.scaled$DistanceToClusterMean <- sapply(1:nrow(qc.scaled),  function(x) {
    this.knns <- knns$nn.index[x,];
    center <- colMeans(qc.scaled[this.knns, 2:(ncol(qc.scaled)-1)]);
    return(as.numeric(dist(rbind(as.numeric(qc.scaled[x, 2:(ncol(qc.scaled)-1)]), as.numeric(center)))))
})


qc.scaled <- as.data.table(qc.scaled)
qc.scaled$label <- qc$V1
qc.scaled[qc.scaled$DistanceToClusterMean < 0.75,]$label <- ""
e = ecdf(qc.scaled$DistanceToClusterMean)
qc.scaled[,`ecdf(NN_distance)`:=e(DistanceToClusterMean)]

# Plot distance distribution

svglite("knn-distance-distribution.svg", width = 10, height = 10)
ggplot(qc.scaled, aes(DistanceToClusterMean)) +
  stat_ecdf() + theme_bw() + geom_text_repel(data = qc.scaled, aes(DistanceToClusterMean,`ecdf(NN_distance)`, label=label ),max.overlaps = 50,size=2)
dev.off()

# write table of metrics.....
mediodDT <- data.table(sample_id=rownames(dd),med_cluster=dd[,"cluster"])
knndistDT <- data.table(sample_id=qc.scaled$V1,mean_dist=qc.scaled$DistanceToClusterMean,ecdf=qc.scaled[,`ecdf(NN_distance)`],n_ref=qc.scaled$n_ref)
fwrite(x=mediodDT[knndistDT,on=.(sample_id)],
       file="knn-kmedoid-stats.txt",col.names=T, row.names=F,quote=F,sep="\t")
