# R script
# David T Severson
# Scope: Functions in BEARscc used to evaluate clusters labels in \\ the context of the cluster labels identified in noise perturbation clusters.

###underlying function calculating cluster stability###############

##create cell association matrices for compute consensus
create_cm<-function(z, names){
  names(z)<-names
  sapply(z, function(x){sapply(z, function(y){as.numeric(x==y)})})
}

##functions for cluster stability from estimate_noise.R
###x is the consensus matrix
###y is a dataframe of cluster identities for various k
calculate_cluster_metrics_by_cluster<-function(i,y,x){
  in.cluster<-as.character(names(y[y==i]))
  out.cluster<-as.character(names(y[!(y==i)]))
  if (length(in.cluster)>1){
    inside.stability<-(sum(x[in.cluster,in.cluster, drop=FALSE])-length(in.cluster))/(length(in.cluster)^2-length(in.cluster))
  } else {
    inside.stability<-1
  }
  promiscuity.matrix<-as.matrix(x[in.cluster,out.cluster, drop=FALSE])
  if (length(in.cluster) >= length(out.cluster)){
    outside.promiscuity<-mean(promiscuity.matrix)
  } else {
    prom.mean.dt<-data.table(data.frame(colMeans(promiscuity.matrix)), keep.rownames = TRUE)
    colnames(prom.mean.dt)[2]<-"mean.prom"
    prom.mean.dt<-prom.mean.dt[order(-mean.prom)]
    prom.mean.dt<-prom.mean.dt[1:length(in.cluster),]
    outside.promiscuity<-mean(promiscuity.matrix[,as.character(prom.mean.dt$rn), drop=FALSE])
  }
  score<-inside.stability-outside.promiscuity
  metrics.dt<-data.table(data.frame(score=score, Promiscuity=outside.promiscuity, Stability=inside.stability, size=length(in.cluster)))
  metrics.dt
}

calculate_cluster_metrics<-function(y,x,z){
  names(y)<-z
  cluster_names<-as.vector(as.character(unique(y)))
  names(cluster_names)<-as.character(unique(y))
  metric_list<-lapply(cluster_names, `calculate_cluster_metrics_by_cluster`, x=x, y=y)
  metric_list<-melt(metric_list, id.vars=c("size"))
  colnames(metric_list)[c(2,4)]<-c("metric","rn")
  metric_list
}

calculate_cell_metrics_by_cluster<-function(i, y, x){
  in.cluster<-as.character(names(y[y==i]))
  out.cluster<-as.character(names(y[!(y==i)]))
  if (length(in.cluster) < length(out.cluster)){
    if (length(in.cluster)==1){
      inside.stability<-1
      outside.promiscuity<-max(x[in.cluster, out.cluster])
      score<-1-outside.promiscuity
      cluster_metric.dt<-data.frame(cell=as.character(in.cluster[1]),score=score, Promiscuity=outside.promiscuity, Stability=inside.stability, cluster.size=length(in.cluster))
      cluster_metric.dt$cluster.size<-length(in.cluster)
    } else {
      stability.dt<-data.table(data.frame(x[in.cluster,in.cluster, drop=FALSE]), keep.rownames = TRUE)
      stability.dt<-stability.dt[,.(Stability = (rowSums(.SD)-1)/(length(in.cluster)-1)), keyby=rn]
      promiscuity.dt<-data.table(data.frame(x[in.cluster,out.cluster, drop=FALSE]), keep.rownames = TRUE)
      promiscuity.dt<-promiscuity.dt[,.(Promiscuity = apply(.SD,1, function(x, y=length(in.cluster)) {mean(sort(x, decreasing = TRUE)[1:y])})), keyby=rn]
      cluster_metric.dt<-stability.dt[promiscuity.dt]
      colnames(cluster_metric.dt)[1]<-"cell"
      cluster_metric.dt$score<-cluster_metric.dt$Stability-cluster_metric.dt$Promiscuity
      cluster_metric.dt$cluster.size<-length(in.cluster)
    }
  } else if (length(in.cluster) >= length(out.cluster)) {
    if (length(in.cluster)==dim(x)[1]) {
      stability.dt<-data.table(data.frame(x[in.cluster,in.cluster, drop=FALSE]), keep.rownames = TRUE)
      stability.dt<-stability.dt[,.(Stability = (rowSums(.SD)-1)/(length(in.cluster)-1)), by=rn]
      stability.dt$Promiscuity<-0
      colnames(stability.dt)[1]<-"cell"
      cluster_metric.dt<-stability.dt
      cluster_metric.dt$cluster.size<-length(in.cluster)
      cluster_metric.dt$score<-cluster_metric.dt$Stability-cluster_metric.dt$Promiscuity
    } else {
      stability.dt<-data.table(data.frame(x[in.cluster,in.cluster,drop=FALSE]), keep.rownames = TRUE)
      stability.dt<-stability.dt[,.(Stability = (rowSums(.SD)-1)/(length(in.cluster)-1)), keyby=rn]
      promiscuity.dt<-data.table(data.frame(x[in.cluster,out.cluster, drop=FALSE]), keep.rownames = TRUE)
      promiscuity.dt<-promiscuity.dt[,.(Promiscuity = rowMeans(.SD)), keyby=rn]
      cluster_metric.dt<-stability.dt[promiscuity.dt]
      cluster_metric.dt$score<-cluster_metric.dt$Stability-cluster_metric.dt$Promiscuity
      cluster_metric.dt$cluster.size<-length(in.cluster)
      colnames(cluster_metric.dt)<-c("cell", "Stability", "Promiscuity", "score", "cluster.size")
    }
  }
  cluster_metric.dt
}

calculate_cell_metrics<-function(y,x,z){
  names(y)<-z
  cluster_names<-as.vector(as.character(unique(y)))
  names(cluster_names)<-as.character(unique(y))
  cell_metric_list<-data.table(melt(lapply(cluster_names, `calculate_cell_metrics_by_cluster`, x=x, y=y), id.vars=c("cell", "cluster.size")))
  colnames(cell_metric_list)[c(3,5)]<-c("metric","rn")
  cell_metric_list
}

###########USER FUNCTIONS##############################
##compute consensus from dataframe of perturbaton reclustering labels
compute_consensus<-function(x, plot=FALSE, file="Rplot", orig.labels=NA, num.clust=5){
  print("Calculating consensus matrix...")
  cm<-lapply(x, `create_cm`,names=rownames(x))
  cm<-Reduce("+",cm)/dim(x)[2]
  print("Calculating noise consensus clusterings...")
  rclust<-hclust(as.dist(abs(1-cm)), method="average")
  cclust<-hclust(as.dist(abs(1-t(cm))), method="average")
  clusts <- cutree(cclust, k = num.clust)
  clusts<-clusts[colnames(cm)]
  nccolors<-c("darkred", "dodgerblue", "blue3", "limegreen", "darkorange4", "darkcyan", "purple", "gold3","#dad1bf", "#5e4400", "#00abff", "#73a3df", "#d6d1c8","#514627", "maroon3", "#87d6ff", "#004eae", "#a4a097", "#efcf72", "#93d6ff", "#304865", "#76726a", "#b99e4e", "#0077c6", "#7da3d6", "#4a463f", "#847035", "#0077cf", "#55749c", "#e0d0ac", "#52461f", "#0f76be", "#84a2cd", "#ddd0b5", "#f2ce5c", "#004a7d", "#0050bf", "#ab9f86", "#ba9e43", "#9dd5ff", "#5d7393", "#a8a08f", "#85702c", "#31a6ff", "#38485e", "#797161", "#544516", "#004a85", "#b6d4ff", "#e6d099", "#f5ce44", "#0078d8", "#0052d0", "#e3d0a3", "#bd9e2b", "#0079e1", "#8ba2c4", "#ae9f7d", "#867020", "#004b8d", "#63738b", "#7c7159", "#8a7000", "#154975", "#3e4756", "#4d4637", "#fbce00", "#a6d5ff", "#bcd3f7", "#ebcf86", "#bf9d00", "#007ae9", "#c2d3ed", "#e9cf8f", "#564500", "#2f76b5", "#92a2bb", "#b39f6a", "#d1d1d1", "#004c95", "#697383", "#b19f74", "#a0a0a0", "#007cff", "#c8d2e4", "#807148", "#727272", "#5ba5f1", "#0054e1", "#7e7150", "#474747", "#4075ad", "#97a1b2", "#4f462f", "#00ddff", "#004c9d", "#9ca1a9", "#edcf7c", "#29d9ff", "#26496d", "#6e727a", "#b79e58", "#50d8ff", "#aed4ff", "#43474e", "#b59e61", "#67d7ff", "#69a4e8", "#cdd2db", "#82703f", "#79d7ff", "#4c74a4")
  newclusters<-as.character(clusts[colnames(cm)])
  if (plot==TRUE && is.na(orig.labels)){
    print("Creating plot without original cluster labels, because none were provided.")
    aheatmap(cm,treeheight=80, annCol = list("Consensus Clusters"=newclusters),labCol=NA, labRow=" ", Rowv = rclust, Colv=cclust)
  } else if (plot==TRUE && !(is.na(orig.labels))) {
    print("Creating plot with original cluster labels, because they were provided")
    oldclusters<-as.character(orig.labels[colnames(cm),])
    aheatmap(cm,treeheight=80, annColors =list("Consensus Clusters"=rev(nccolors), "Original Clusters"=nccolors) ,annCol = list("Consensus Clusters"=newclusters, "Original Clusters"=oldclusters), labCol=NA, labRow=" ", Rowv = rclust, Colv=cclust)
  }
  clusts<-data.table(data.frame(clusts), keep.rownames = TRUE)[order(clusts)]
  list(consensus_matrix=cm,clusts=clusts)
}

cluster_consensus<-function(x,y, method="complete"){
  cclust<-hclust(as.dist(abs(1-t(x))), method=method)
  clusts<-data.frame(cutree(cclust, k=y))
  clusts
}

##seperately plot cell and cluster noise robustness index
# y is a data frame of the final clustering results perturbed labels
# x is a data frame of the consensus matrix. need to make robust to bug when names start with numbers...
report_cell_metrics<-function(y,x){
  z<-rownames(x)
  y$X1<-rep("1",dim(y)[1])
  y<-data.table(y)
  cell.scores<-list()
  cell.scores<-lapply(y, `calculate_cell_metrics`,x=x, z=z)
  cell.scores<-data.table(melt(cell.scores, id.var=c("rn","cell","cluster.size", "metric")))
  cell.scores$L1<-as.factor(gsub("^X","", cell.scores$L1))
  cell.scores<-cell.scores[,-5, with=FALSE]
  cell.scores
}

###x is the consensus matrix with rows and columns as cell names
###y is a dataframe of cluster identities for various k in columns and rownames as cells
report_cluster_metrics<-function(y,x, weighted_mean=FALSE,plot=FALSE, file="Rplot"){
  print(weighted_mean)
  z<-rownames(y)
  y<-data.table(y)
  cluster.scores<-list()
  oneclust.stability<-(sum(x)-dim(x)[1])/dim(x)[1]^2
  cluster.scores$X1<-melt(data.table(rn="1",score=oneclust.stability,Promiscuity=0,Stability=oneclust.stability, size=dim(x)[1]), id.vars=c("rn", "size"))
  colnames(cluster.scores$X1)[3]<-"metric"
  cluster.scores<-append(cluster.scores,lapply(y, `calculate_cluster_metrics`,x=x, z=z))
  cluster.scores.m<-data.table(melt(cluster.scores, id.var=c("rn","size","metric")))
  cluster.scores.m<-cluster.scores.m[,-4, with=FALSE]
  cluster.scores.m$L1<-gsub("^X","", cluster.scores.m$L1)
  cluster.scores.m$singlet<-cluster.scores.m$size==1
  cluster.scores.m$singlet<-factor(gsub(FALSE,"Cell number > 1", gsub(TRUE, "Cell number = 1", cluster.scores.m$singlet)),levels = c("Cell number > 1","Cell number = 1"))
  if (weighted_mean==TRUE){
    mean.overall<-cluster.scores.m[,sum(value*size/as.numeric(cluster.scores.m[L1=="1",]$size[1])), by=c("L1", "metric")]
  }else {
    mean.overall<-cluster.scores.m[,mean(value), by=c("L1", "metric")]
  }
  colnames(mean.overall)[3]<-"Overall.mean"
  setkeyv(cluster.scores.m,c("L1", "metric"))
  setkeyv(mean.overall, c("L1", "metric"))
  cluster.scores.m<-cluster.scores.m[mean.overall]
  numbers<-as.character(sort(unique(na.omit(as.numeric(gsub("^X","",cluster.scores.m$L1))))))
  other<-unique(cluster.scores.m$L1[!(gsub("X","",cluster.scores.m$L1)%in%numbers)])
  level<-as.character(c(numbers,other))
  cluster.scores.m$L1<-factor(cluster.scores.m$L1, levels = level)
  if (plot==TRUE){
    g<-ggplot(data=cluster.scores.m, aes(x=L1, y=value))+geom_boxplot(color="grey20", fill="grey85")+
      guides(fill="none", color="none")+xlab("\nCluster number")+ylab("Cluster metric\n")+facet_wrap(~metric, ncol=1,nrow=3, scales = "free")+
      ggtitle("Cluster technical consensus scores for various clustering\n")+labs(alpha="Cluster size")+
      geom_point(aes(x=factor(L1), y=Overall.mean),fill="red", size=2.5, color="black", pch=23)+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
    savefile<-paste(file, "_cluster_scores.pdf", sep="")
    ggsave(savefile, plot=g, device="pdf", width=8.5, height=11, units = "in")
  }
  cluster.scores.m
}
