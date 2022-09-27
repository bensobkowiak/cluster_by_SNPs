


cluster_by_SNPs<-function(SNPdistance, thresholds=c(5,10,20,50,100),output){
  options(stringsAsFactors = F)
  distMatrix<-read.table(SNPdistance,check.names = F)
  distMatrix[lower.tri(distMatrix,diag = T)] = NA
  distmat<-reshape2::melt(as.matrix(distMatrix), varnames = c('row', 'col'), na.rm = TRUE)
  row.names(distmat)<-NULL
  thresholds<-thresholds
  results<-data.frame(matrix(NA,ncol = length(thresholds)+1,nrow=nrow(distMatrix)))
  results[,1]<-row.names(distMatrix)
  
  for (j in 1:length(thresholds)){
    close_relate<-distmat[which(as.numeric(distmat[,3])<=thresholds[j]),]
    names<-unique(c(as.character(close_relate[,1]),as.character(close_relate[,2])))
    cluster_results<-as.data.frame(matrix(0,ncol = 2, nrow = length(names)))
    cluster_results[,1]<-as.character(names)
    cluster_results[,2]<-1:length(names)
    for (i in 1:nrow(cluster_results)){
      newclose<-close_relate[which(close_relate$row==cluster_results[i,1] | close_relate$col==cluster_results[i,1]),]
      clusters<-unique(cluster_results[which(cluster_results[,1] %in% 
                                        unique(c(as.character(newclose[,1]),as.character(newclose[,2])))),2])
      clusters<-clusters[!is.na(clusters)]
      cluster_results[which(cluster_results[,1] %in% 
                              unique(c(as.character(newclose[,1]),as.character(newclose[,2])))),2]<-cluster_results[i,2]
      cluster_results[which(cluster_results[,2] %in% clusters),2]<-cluster_results[i,2]
      close_relate<-close_relate[-which(close_relate$row==cluster_results[i,1] | close_relate$col==cluster_results[i,1]),]
    }
    row.names(cluster_results)<-NULL
    for (k in 1:nrow(cluster_results)){
      results[which(cluster_results[k,1]==results[,1]),j+1]<-cluster_results[k,2]
    }
  }
  colnames(results)<-c("SampleID",paste0("SNPs",thresholds))
  results<-results[order(results[,2]),]
  
  for (i in 2:ncol(results)){ ## rename clusters
    clustnames<-unique(results[,i])
    for (j in 1:length(clustnames)){
      if (!is.na(clustnames[j])){
        results[which(results[,i]==clustnames[j]),i]<-paste0("cluster_",j)
      }
    }
  }
  write.csv(results,paste0(output,".csv"),row.names = F)
}
