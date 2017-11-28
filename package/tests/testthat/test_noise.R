#R
#Author: David T Severson
#Test functions for BEARscc package

## load example data and run BEARscc workflow for tests ##
cat("\n")
print("Running BEARscc on example data to test for correct installation.")
cat("\n")
data("BEARscc_examples")
results<-estimate_noiseparameters(ERCC.counts.df, data.counts.df,ERCC.meta.df, granularity=30, write.noise.model=FALSE, plot=FALSE, alpha_granularity = 0.25)
test_that("estimation output?", {
  expect_true(is.list(results))
  expect_equal(dim(results$bayes_parameters)[2], dim(results$original.counts)[1]+4)
  expect_match(colnames(results$bayes_parameters)[1], "rn")
  expect_match(colnames(results$bayes_parameters)[2], "k")
  expect_match(colnames(results$bayes_parameters)[3], "Py0givenk")
  expect_match(colnames(results$bayes_parameters)[4], "Counts")
  expect_equal(dim(results$ERCC_parameters)[2],6)
})
cat("\n")
print("Injecting noise for testing...")
noisy_counts.list<-create_noiseinjected_counts(results,n=3)
test_that("noise perturbation output?", {
  test_noisy_counts<-function(x){
    noisy_counts.dt<-x
    expect_true(is.data.table(noisy_counts.dt))
    expect_equal(dim(noisy_counts.dt)[2], dim(results$original.counts)[2]+1)
    expect_equal(dim(noisy_counts.dt)[1], dim(results$original.counts)[1])
    expect_match(colnames(noisy_counts.dt)[1], "GENE_ID")
  }
  lapply(noisy_counts.list, `test_noisy_counts`)
})
cat("\n")
print("Creating consensus matrix for testing...")
recluster<-function(x){
  x<-data.frame(x, row.names = "GENE_ID")
  scramble<-sample(colnames(x), size=length(colnames(x)), replace=FALSE)
  x<-x[,scramble]
  clust<-hclust(dist(t(x),method="euclidean"),method="complete")
  clust<-cutree(clust,2)
  data.frame(clust)
}
clusters.list<-lapply(noisy_counts.list, `recluster`)
clusters.df<-do.call("cbind", clusters.list)
colnames(clusters.df)<-names(clusters.list)
#write.table(clusters.df, file="example_clusters.tsv", col.names = TRUE, row.names = TRUE,quote=FALSE, sep = "\t")
noise_consensus<-compute_consensus(clusters.df)
test_that("consensus matrix result?", {
  expect_true(sum(noise_consensus>1)==0)
  expect_true(sum(noise_consensus<0)==0)
  expect_equal(dim(noise_consensus)[1], dim(results$original.counts)[2])
  expect_equal(dim(noise_consensus)[2], dim(results$original.counts)[2])
})
cat("\n")
print("Computing cluster and cell metrics for testing...")
vector<-seq(from=2, to=5, by=1)
BEARscc_clusts.df<-cluster_consensus(noise_consensus,vector)
BEARscc_clusts.df<-cbind(BEARscc_clusts.df, Original=clusters.df$Original_counts)
cluster_scores.dt<-report_cluster_metrics(BEARscc_clusts.df,noise_consensus, plot=FALSE)
cell_scores.dt<-report_cell_metrics(BEARscc_clusts.df, noise_consensus)

test_that("cluster metrics are computing well? ", {
  expect_equal(dim(BEARscc_clusts.df)[1], dim(noise_consensus)[1])
  expect_equal(dim(BEARscc_clusts.df)[2], eval(length(vector)+1))
  expect_equal(dim(cluster_scores.dt)[2], 7)
  expect_true(sum(as.character(unique(cluster_scores.dt$Clustering)) %in% gsub("^X","",colnames(BEARscc_clusts.df)))==5)
  expect_equal(dim(cell_scores.dt)[2], 6)
  expect_true(sum(as.character(unique(cell_scores.dt$Clustering)) %in% gsub("^X","",colnames(BEARscc_clusts.df)))==5)
  expect_true(sum(unique(cell_scores.dt$Cell) %in% rownames(BEARscc_clusts.df))==50)
})

