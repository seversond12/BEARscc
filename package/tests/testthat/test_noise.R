#R
#Author: David T Severson
#Test functions for BEARscc package

## load example data and run BEARscc for tests ##
data.counts.dt<- fread("brain_control_example.tsv")
ERCC.meta.dt<- fread("ERCC.meta.tsv")
ERCC.counts.df<-data.frame(data.counts.dt[GENE_ID%like%"ERCC-",], row.names="GENE_ID")
data.counts.df<-data.frame(data.counts.dt, row.names = "GENE_ID")
ERCC.meta.df<-data.frame(ERCC.meta.dt, row.names="ERCC_ID")
results<-estimate_noiseparameters(ERCC.counts.df,ERCC.meta.df, data.counts.df, granularity=30, write.noise.model=FALSE, plot=FALSE, alpha_granularity = 0.25)
noisy_counts.list<-create_noiseinjected_counts(results)

test_that("estimation output?", {
  expect_true(is.list(results))
  expect_equal(dim(results$bayes_parameters)[2], dim(results$original.counts)[1]+4)
  expect_match(colnames(results$bayes_parameters)[1], "rn")
  expect_match(colnames(results$bayes_parameters)[2], "k")
  expect_match(colnames(results$bayes_parameters)[3], "Py0givenk")
  expect_match(colnames(results$bayes_parameters)[4], "Counts")
  expect_equal(dim(results$ERCC_parameters)[2],6)
})


test_that("noise perturbation output?", {
  test_noisy_counts<-function(x){
    noisy_counts.dt<-x
    expect_true(is.data.table(noisy_counts.dt))
    expect_equal(dim(noisy_counts.dt)[2], dim(results$original.counts)[2]+1)
    expect_equal(dim(noisy_counts.dt)[1], dim(results$original.counts)[1])
    expect_match(colnames(noisy_counts.dt)[1], "GENE_ID")
    expect_true(identical(data.frame(noisy_counts.dt, row.names="GENE_ID"), round(data.frame(noisy_counts.dt, row.names="GENE_ID"))))
  }
  lapply(noisy_counts.list, `test_noisy_counts`)
})


