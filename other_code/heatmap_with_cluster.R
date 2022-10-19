args <- commandArgs(TRUE)
pdf(args[1],width=10,height=10)
data<-read.table(args[2])
data<- as.matrix(data)
library(gplots)
heatmap.2(data,scale="row",key=T,keysize=0.8,trace="none",cexCol=0.5,cexRow=0.5,srtCol=90,offsetRow=0,margins=c(5,10))#col=greenred
#heatmap.2(data,scale="row",key=T,keysize=0.8,trace="none",cexCol=0.5,labRow=NA,srtCol=45,offsetRow=0,margins=c(5,8))#col=greenred
#heatmap.2(data,scale="row",key=T,keysize=1.5,trace="none",cexCol=1,cexRow=0.5,offsetRow=0,margins=c(5,8),Rowv=NA,Colv=NA,dendrogram=c("none"))
#standardization the data, in case of some abnormal data (scale=("none","row","column"))
#if you don't have legend, you can set key=F. keysize is the size of legend
#trace can move line in the heatmap
#cexCol and cexRow is represnet the size of the xlab and ylab
#dendrogram is cluster
