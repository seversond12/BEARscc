# R script
# David T Severson
# Scope: Functions in BEARscc used to evaluate clusters labels in
    # the context of the cluster labels identified in
    #noise perturbation clusters.

###underlying function calculating cluster stability###############

#creates a consensus matrix from a single list of cluster labels
create_cm<-function(cluster_labels, names){
    names(cluster_labels)<-names
    sapply(cluster_labels, function(x){sapply(cluster_labels,
        function(y){as.numeric(x==y)})})
}

calculate_cluster_metrics_by_cluster<-function(clusters, cluster_labels,
consensus_matrix){
    mean.prom<-NULL
    in.cluster<-as.character(names(cluster_labels[cluster_labels==clusters]))
    out.cluster<-as.character(names(
        cluster_labels[!(cluster_labels==clusters)]))
    if (length(in.cluster)>1){
        inside.stability<-(sum(consensus_matrix[in.cluster,in.cluster,
            drop=FALSE])-length(in.cluster))/
            (length(in.cluster)^2-length(in.cluster))
    } else {
        inside.stability<-1
    }
    promiscuity.matrix<-as.matrix(consensus_matrix[in.cluster,out.cluster,
        drop=FALSE])
    if (length(in.cluster) >= length(out.cluster)){
        outside.promiscuity<-mean(promiscuity.matrix)
    } else {
        prom.mean.dt<-data.table(data.frame(colMeans(promiscuity.matrix)),
            keep.rownames = TRUE)
        colnames(prom.mean.dt)[2]<-"mean.prom"
        prom.mean.dt<-prom.mean.dt[order(-mean.prom)]
        prom.mean.dt<-prom.mean.dt[1:length(in.cluster),]
        outside.promiscuity<-mean(promiscuity.matrix[,
            as.character(prom.mean.dt$rn), drop=FALSE])
    }
    score<-inside.stability-outside.promiscuity
    metrics.dt<-data.table(data.frame(Score=score,
        Promiscuity=outside.promiscuity, Stability=inside.stability,
        size=length(in.cluster)))
    metrics.dt
}

calculate_cluster_metrics<-function(cluster_labels, consensus_matrix,
cluster_names){
    names(cluster_labels)<-cluster_names
    clusters<-as.vector(as.character(unique(cluster_labels)))
    names(clusters)<-as.character(unique(cluster_labels))
    metric_list<-lapply(clusters, `calculate_cluster_metrics_by_cluster`,
        consensus_matrix=consensus_matrix, cluster_labels=cluster_labels)
    metric_list<-melt(metric_list, id.vars=c("size"))
    colnames(metric_list)[c(2,4)]<-c("metric","rn")
    metric_list
}

calculate_cell_metrics_by_cluster<-function(cluster_names, cluster_labels,
consensus_matrix){
    .<-rn<-NULL
    in.cluster<-as.character(names(
        cluster_labels[cluster_labels==cluster_names]))
    out.cluster<-as.character(names(
        cluster_labels[!(cluster_labels==cluster_names)]))
    if (length(in.cluster) < length(out.cluster)){
        if (length(in.cluster)==1){
            inside.stability<-1
            outside.promiscuity<-max(consensus_matrix[in.cluster, out.cluster])
            score<-1-outside.promiscuity
            cell_metric.dt<-data.frame(cell=as.character(in.cluster[1]),
            Score=score, Promiscuity=outside.promiscuity,
            Stability=inside.stability, cluster.size=length(in.cluster))
            cell_metric.dt$cluster.size<-length(in.cluster)
        } else {
            stability.dt<-data.table(data.frame(consensus_matrix[in.cluster,
            in.cluster, drop=FALSE]), keep.rownames = TRUE)
            stability.dt<-stability.dt[,
            .(Stability = (rowSums(.SD)-1)/(length(in.cluster)-1)), keyby=rn]
            promiscuity.dt<-data.table(data.frame(
            consensus_matrix[in.cluster,out.cluster, drop=FALSE]),
            keep.rownames = TRUE)
            promiscuity.dt<-promiscuity.dt[,.(Promiscuity = apply(.SD, 1,
            function(x, y=length(in.cluster)) {mean(sort(x,
            decreasing = TRUE)[1:y])})), keyby=rn]
            cell_metric.dt<-stability.dt[promiscuity.dt]
            colnames(cell_metric.dt)[1]<-"cell"
            cell_metric.dt$Score<-cell_metric.dt$Stability-
            cell_metric.dt$Promiscuity
            cell_metric.dt$cluster.size<-length(in.cluster)
        }
    } else if (length(in.cluster) >= length(out.cluster)) {
        if (length(in.cluster)==dim(consensus_matrix)[1]) {
            stability.dt<-data.table(data.frame(consensus_matrix[in.cluster,
                in.cluster, drop=FALSE]), keep.rownames = TRUE)
            stability.dt<-stability.dt[,
                .(Stability = (rowSums(.SD)-1)/(length(in.cluster)-1)), by=rn]
            stability.dt$Promiscuity<-0
            colnames(stability.dt)[1]<-"cell"
            cell_metric.dt<-stability.dt
            cell_metric.dt$cluster.size<-length(in.cluster)
            cell_metric.dt$Score<-cell_metric.dt$Stability-
            cell_metric.dt$Promiscuity
        } else {
            stability.dt<-data.table(data.frame(consensus_matrix[in.cluster,
                in.cluster,drop=FALSE]), keep.rownames = TRUE)
            stability.dt<-stability.dt[,
                .(Stability = (rowSums(.SD)-1)/(length(in.cluster)-1)), keyby=rn]
            promiscuity.dt<-data.table(data.frame(consensus_matrix[in.cluster,
                out.cluster, drop=FALSE]), keep.rownames = TRUE)
            promiscuity.dt<-promiscuity.dt[,.(Promiscuity = rowMeans(.SD)),
                keyby=rn]
            cell_metric.dt<-stability.dt[promiscuity.dt]
                cell_metric.dt$Score<-cell_metric.dt$Stability-
                cell_metric.dt$Promiscuity
            cell_metric.dt$cluster.size<-length(in.cluster)
            colnames(cell_metric.dt)<-c("cell", "Stability",
                "Promiscuity", "Score", "cluster.size")
        }
    }
    cell_metric.dt
}

calculate_cell_metrics<-function(cluster_labels,consensus_matrix,clusters){
    names(cluster_labels)<-clusters
    cluster_names<-as.vector(as.character(unique(cluster_labels)))
    cluster_names<-as.character(unique(cluster_labels))
    cell_metric_list<-data.table(melt(lapply(cluster_names,
        `calculate_cell_metrics_by_cluster`, consensus_matrix=consensus_matrix,
        cluster_labels=cluster_labels), id.vars=c("cell", "cluster.size")))
    colnames(cell_metric_list)[c(3,5)]<-c("metric","rn")
    cell_metric_list
}


plot_cluster_metrics<-function(cluster.scores.m, file){
    L1<-value<-Overall.mean<-NULL
    g<-ggplot(data=cluster.scores.m, aes(x=L1, y=value))+
        geom_boxplot(color="grey20", fill="grey85")+
        guides(fill="none", color="none")+xlab("\nCluster number")+
        ylab("Cluster metric\n")+
        facet_wrap(~metric, ncol=1,nrow=3, scales = "free")+
        ggtitle("Cluster technical consensus scores
        for various clustering\n")+labs(alpha="Cluster size")+
        geom_point(aes(x=factor(L1), y=Overall.mean),fill="red", size=2.5,
        color="black", pch=23)+
        annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
        annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    savefile<-paste(file, "_cluster_scores.pdf", sep="")
    ggsave(savefile, plot=g, device="pdf", width=8.5,
        height=11, units = "in")
}

###########USER FUNCTIONS##############################
#user function for computing a consensus matrix from cluster labels
compute_consensus<-function(cluster_labels){
    orig.clusts<-cluster_labels[,length(colnames(cluster_labels))]
    cluster_labels<-cluster_labels[,-length(colnames(cluster_labels))]
    cm<-lapply(cluster_labels, `create_cm`,names=rownames(cluster_labels))
    cm<-Reduce("+",cm)/dim(cluster_labels)[2]
    cm
}

#user function for clustering the noise consensus matrix
cluster_consensus<-function(consensus_matrix, cluster_num, method="complete"){
    cclust<-hclust(as.dist(abs(1-t(consensus_matrix))), method=method)
    clusts<-data.frame(cutree(cclust, k=cluster_num))
    clusts
}

##user function for computing BEARscc cell metrics
report_cell_metrics<-function(cluster_labels, consensus_matrix){
    clusters<-rownames(cluster_labels)
    cluster_labels$X1<-rep("1",dim(cluster_labels)[1])
    cluster_labels<-data.table(cluster_labels)
    cell.scores<-list()
    cell.scores<-lapply(cluster_labels, `calculate_cell_metrics`,
        consensus_matrix=consensus_matrix, clusters=clusters)
    cell.scores<-data.table(melt(cell.scores, id.var=c("rn", "cell",
        "cluster.size", "metric")))
    cell.scores$L1<-as.factor(gsub("^X","", cell.scores$L1))
    cell.scores<-cell.scores[,-5, with=FALSE]
    colnames(cell.scores)<-c("Cluster identity", "Cell",
        "Cluster size", "Metric", "Value", "Clustering")
    data.frame(cell.scores)
}

##user function for computing BEARscc cluster metrics
report_cluster_metrics<-function(cluster_labels, consensus_matrix,
weighted_mean=FALSE, plot=FALSE, file="Rplot"){
    values<-value<-size<-L1<-Overall.mean<-NULL
    cluster_names<-rownames(cluster_labels)
    cluster_labels<-data.table(cluster_labels)
    cluster.scores<-list()
    oneclust.stability<-(sum(consensus_matrix)-dim(consensus_matrix)[1])/
        dim(consensus_matrix)[1]^2
    cluster.scores$X1<-melt(data.table(rn="1", Score=oneclust.stability,
        Promiscuity=0, Stability=oneclust.stability,
        size=dim(consensus_matrix)[1]), id.vars=c("rn", "size"))
    colnames(cluster.scores$X1)[3]<-"metric"
    cluster.scores<-append(cluster.scores, lapply(cluster_labels,
        `calculate_cluster_metrics`, consensus_matrix=consensus_matrix,
        cluster_names=cluster_names))
    cluster.scores.m<-data.table(melt(cluster.scores,
        id.var=c("rn","size","metric")))
    cluster.scores.m<-cluster.scores.m[,-4, with=FALSE]
    cluster.scores.m$L1<-gsub("^X","", cluster.scores.m$L1)
    cluster.scores.m$singlet<-cluster.scores.m$size==1
    cluster.scores.m$singlet<-factor(gsub(FALSE,"Cell number > 1",
        gsub(TRUE, "Cell number = 1", cluster.scores.m$singlet)),
        levels = c("Cell number > 1","Cell number = 1"))
    if (weighted_mean==TRUE){
        mean.overall<-cluster.scores.m[,
            sum(value*size/as.numeric(cluster.scores.m[L1=="1",]$size[1])),
            by=c("L1", "metric")]
    }else {
        mean.overall<-cluster.scores.m[,mean(value), by=c("L1", "metric")]
    }
    colnames(mean.overall)[3]<-"Overall.mean"
    setkeyv(cluster.scores.m,c("L1", "metric"))
    setkeyv(mean.overall, c("L1", "metric"))
    cluster.scores.m<-cluster.scores.m[mean.overall]
    numbers<-as.character(sort(unique(suppressWarnings(
        na.omit(as.numeric(gsub("^X","",cluster.scores.m$L1)))))))
    other<-unique(cluster.scores.m$L1[!(gsub("X", "",
        cluster.scores.m$L1)%in%numbers)])
    level<-as.character(c(numbers,other))
    cluster.scores.m$L1<-factor(cluster.scores.m$L1, levels = level)
    if (plot==TRUE){
        plot_cluster_metrics(cluster.scores.m=cluster.scores.m, file=file)
    }
    colnames(cluster.scores.m)<-c("Cluster identity", "Cluster size",
        "Metric", "Value", "Clustering", "Singlet", "Clustering Mean")
    data.frame(cluster.scores.m)
}
