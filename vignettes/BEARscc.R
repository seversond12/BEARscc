## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----prepare_libaries, eval=TRUE, include=FALSE----------------------------
library("data.table")
library("SingleCellExperiment")

## ----install_BEARscc, include=TRUE, eval=FALSE-----------------------------
#  install.packages('../inst/BEARscc_0.99.3.tar.gz', repos = NULL, type="source")

## ----load_data, eval=TRUE--------------------------------------------------
library("BEARscc")
data("BEARscc_examples")

## ----display_data, eval=TRUE-----------------------------------------------
head(ERCC.counts.df[,1:2])

head(data.counts.df[,1:2]) 

head(ERCC.meta.df)

BEAR_examples.sce 

## ----create_SCEList, eval=TRUE---------------------------------------------
BEAR.se <- SummarizedExperiment(list(counts=as.matrix(data.counts.df)))
BEAR_examples.sce<-as(BEAR.se, "SingleCellExperiment")
metadata(BEAR_examples.sce)<-list(spikeConcentrations=ERCC.meta.df)
assay(BEAR_examples.sce, "observed_expression")<-counts(BEAR_examples.sce)
isSpike(BEAR_examples.sce, "ERCC_spikes")<-grepl("^ERCC-", 
    rownames(BEAR_examples.sce))

## ----estimate_noise, eval=TRUE---------------------------------------------
BEAR_examples.sce <- estimate_noiseparameters(BEAR_examples.sce,
    max_cumprob=0.9999, alpha_resolution = 0.1, bins=10,
    write.noise.model=FALSE, file="BEAR_examples")

## ----simulate_replicates, eval=TRUE----------------------------------------
    BEAR_examples.sce <- simulate_replicates(BEAR_examples.sce,
        max_cumprob=0.9999, n=10)

## ----estimate_noise_HPC, eval=FALSE, include=TRUE--------------------------
#  BEAR_examples.sce <- estimate_noiseparameters(BEAR_examples.sce,
#                                      write.noise.model=TRUE,
#                                      file="tutorial_example",
#                                      model_view=c("Observed","Optimized"))

## ----HPC_example, eval=FALSE, include=TRUE---------------------------------
#  library("BEARscc")
#  
#  #### Load data ####
#  ITERATION<-commandArgs(trailingOnly=TRUE)[1]
#  counts.df<-read.delim("tutorial_example_counts4clusterperturbation.xls")
#  #filter out zero counts to speed up algorithm
#  counts.df<-counts.df[rowSums(counts.df)>0,]
#  probs4detection<-fread("tutorial_example_bayesianestimates.xls")
#  parameters<-fread("tutorial_example_parameters4randomize.xls")
#  
#  #### Simulate replicates ####
#  counts.error<-HPC_simulate_replicates(counts_matrix=counts.df,
#      dropout_parameters=dropout_parameters,
#      spikein_parameters=spikein_parameters)
#  write.table(counts.error, file=paste("simulated_replicates/",
#      paste(ITERATION,"sim_replicate_counts.txt",sep="_"),
#      sep=""), quote =FALSE, row.names=TRUE)
#  #########

## ----recluster, include=TRUE, eval=TRUE------------------------------------
recluster <- function(x) {
    x <- data.frame(x)
    scramble <- sample(colnames(x), size=length(colnames(x)), replace=FALSE)
    x <- x[,scramble]
    clust <- hclust(dist(t(x), method="euclidean"),method="complete")
    clust <- cutree(clust,2)
    data.frame(clust)
}

## ----original_clusters, include=TRUE, eval=TRUE----------------------------
OG_clusters<-recluster(data.frame(assay(BEAR_examples.sce, 
    "observed_expression")))
colnames(OG_clusters)<-"Original"

## ----recluster_go, include=TRUE, evalue=TRUE-------------------------------
cluster.list<-lapply(metadata(BEAR_examples.sce)$simulated_replicates, 
    `recluster`)
clusters.df<-do.call("cbind", cluster.list)
colnames(clusters.df)<-names(cluster.list)

## ----compute_consensus, eval=TRUE, include=TRUE----------------------------
noise_consensus <- compute_consensus(clusters.df)
head(noise_consensus[,1:3], n=3)

## ----plot_noise_consensus, eval=TRUE, include=TRUE-------------------------
library("NMF")
aheatmap(noise_consensus, breaks=0.5)

## ----cluster_consensus, include=TRUE, eval=TRUE----------------------------
vector <- seq(from=2, to=5, by=1)
BEARscc_clusts.df <- cluster_consensus(noise_consensus,vector)

## ----add_original, include=TRUE, eval=TRUE---------------------------------
BEARscc_clusts.df <- cbind(BEARscc_clusts.df, 
    Original=OG_clusters)

## ----compute_cluster_scores, include=TRUE, eval=TRUE-----------------------
cluster_scores.df <- report_cluster_metrics(BEARscc_clusts.df, 
    noise_consensus, plot=FALSE)
head(cluster_scores.df, n=5)

## ----plot_clusterscores_original, include=TRUE, eval=TRUE------------------
library("ggplot2")
library("cowplot")
original_scores.df<-cluster_scores.df[
    cluster_scores.df$Clustering=="Original",]
ggplot(original_scores.df[original_scores.df$Metric=="Score",], 
    aes(x=`Cluster.identity`, y=Value) )+
    geom_bar(aes(fill=`Cluster.identity`), stat="identity")+
    xlab("Cluster identity")+ylab("Cluster score")+
    ggtitle("Original cluster scores")+guides(fill=FALSE)

## ----plot_clusterstability_original, include=TRUE, eval=TRUE---------------
ggplot(original_scores.df[original_scores.df$Metric=="Stability",], 
    aes(x=`Cluster.identity`, y=Value) )+
    geom_bar(aes(fill=`Cluster.identity`), stat="identity")+
    xlab("Cluster identity")+ylab("Cluster stability")+
    ggtitle("Original cluster stability")+guides(fill=FALSE)

## ----plot_clusterpromiscuity_original, include=TRUE, eval=TRUE-------------
ggplot(original_scores.df[original_scores.df$Metric=="Promiscuity",], 
    aes(x=`Cluster.identity`, y=Value) )+
    geom_bar(aes(fill=`Cluster.identity`), stat="identity")+
    xlab("Cluster identity")+ylab("Cluster promiscuity")+
    ggtitle("Original cluster promiscuity")+guides(fill=FALSE)

## ----compute_cell_scores, include=TRUE, eval=TRUE--------------------------
cell_scores.df <- report_cell_metrics(BEARscc_clusts.df, noise_consensus)
head(cell_scores.df, n=4)

## ----plot_cellscore_original, include=TRUE, eval=TRUE----------------------
original_cell_scores.df<-cell_scores.df[cell_scores.df$Clustering=="Original",]
ggplot(original_cell_scores.df[original_cell_scores.df$Metric=="Score",],
    aes(x=factor(`Cluster.identity`), y=Value) )+
    geom_jitter(aes(color=factor(`Cluster.identity`)), stat="identity")+
    xlab("Cluster identity")+ylab("Cell scores")+
    ggtitle("Original clustering cell scores")+guides(color=FALSE)

## ----choosing_k, include=TRUE, eval=TRUE-----------------------------------
ggplot(cluster_scores.df[cluster_scores.df$Metric=="Score",], 
    aes(x=`Clustering`, y=Value) )+ geom_boxplot(aes(fill=Clustering,  
    color=`Clustering`))+ xlab("Clustering")+
    ylab("Cluster score distribution")+
    ggtitle("Original cluster scores for each clustering")+
    guides(fill=FALSE, color=FALSE)

