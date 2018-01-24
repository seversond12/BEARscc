#R
#Author: David T Severson
#Test functions for BEARscc package

## load example data and run BEARscc workflow for tests ##
library("SingleCellExperiment")
cat("\n")

print("Running BEARscc on example data to test for correct installation.")
cat("\n")

data("BEARscc_examples")
results<-estimate_noiseparameters(BEAR_examples.sce,
    write.noise.model=FALSE, plot=FALSE, alpha_resolution=0.25)
test_that("estimation output?", {
    expect_equal(dim(metadata(results)$dropout_parameters)[2],
        dim(assay(BEAR_examples.sce, "observed_expression"))[1]+3)
    expect_match(colnames(metadata(results)$dropout_parameters)[1],
        "transcripts")
    expect_match(colnames(metadata(results)$dropout_parameters)[2],
        "dropouts_given_transcripts")
    expect_match(colnames(metadata(results)$dropout_parameters)[3],
        "counts")
    expect_equal(dim(metadata(results)$spikein_parameters)[2],6)
})
cat("\n")

print("Injecting noise for testing...")
results<-simulate_replicates(results,n=3)
test_that("noise perturbation output?", {
    test_noisy_counts<-function(x){
        noisy_counts<-x
        expect_true(is.matrix(noisy_counts))
        expect_equal(dim(noisy_counts)[2],
            dim(assay(results, "observed_expression"))[2])
        expect_equal(dim(noisy_counts)[1],
            dim(assay(results, "observed_expression"))[1])
    }
    lapply(metadata(results)$simulated_replicates, `test_noisy_counts`)
})
noisy_counts.list<-metadata(results)$simulated_replicates
cat("\n")

print("Creating consensus matrix for testing...")
recluster<-function(x){
    x<-data.frame(x)
    scramble<-sample(colnames(x), size=length(colnames(x)), replace=FALSE)
    x<-x[,scramble]
    clust<-hclust(dist(t(x),method="euclidean"),method="complete")
    clust<-cutree(clust,2)
    data.frame(clust)
}
clusters.list<-lapply(noisy_counts.list, `recluster`)
clusters.df<-do.call("cbind", clusters.list)
colnames(clusters.df)<-names(clusters.list)
#write.table(clusters.df, file="example_clusters.tsv", col.names = TRUE,
    #row.names = TRUE,quote=FALSE, sep = "\t")
noise_consensus<-compute_consensus(clusters.df)
test_that("consensus matrix result?", {
    expect_true(sum(noise_consensus>1)==0)
    expect_true(sum(noise_consensus<0)==0)
    expect_equal(dim(noise_consensus)[1], dim(assay(results,
        "observed_expression"))[2])
    expect_equal(dim(noise_consensus)[2], dim(assay(results,
        "observed_expression"))[2])
})
cat("\n")

print("Computing cluster and cell metrics for testing...")
OG_clusters<-recluster(data.frame(assay(BEAR_examples.sce,
    "observed_expression")))
vector<-seq(from=2, to=5, by=1)
BEARscc_clusts.df<-cluster_consensus(noise_consensus,vector)
BEARscc_clusts.df<-cbind(BEARscc_clusts.df, Original=OG_clusters)
cluster_scores.dt<-report_cluster_metrics(BEARscc_clusts.df, noise_consensus,
    plot=FALSE)
cell_scores.dt<-report_cell_metrics(BEARscc_clusts.df, noise_consensus)

test_that("cluster metrics are computing well? ", {
    expect_equal(dim(BEARscc_clusts.df)[1], dim(noise_consensus)[1])
    expect_equal(dim(BEARscc_clusts.df)[2], eval(length(vector)+1))
    expect_equal(dim(cluster_scores.dt)[2], 7)
    expect_true(sum(as.character(unique(cluster_scores.dt$Clustering))
        %in% gsub("^X","",colnames(BEARscc_clusts.df)))==5)
    expect_equal(dim(cell_scores.dt)[2], 6)
    expect_true(sum(as.character(unique(cell_scores.dt$Clustering))
        %in% gsub("^X","",colnames(BEARscc_clusts.df)))==5)
    expect_true(sum(unique(cell_scores.dt$Cell) %in% rownames(BEARscc_clusts.df))==50)
})

