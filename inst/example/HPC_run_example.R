#!/usr/bin/env Rscript
#Author: David T Severson
#Scope: Example to run `prepare_probabilities()` on  a high perfromance
    #compute cluster with parallel performance

library("data.table")
library("BEARscc")
library("parallel")

#### Load data ####
ITERATION<-commandArgs(trailingOnly=TRUE)[1]
no_cores<-4
counts.dt<-fread("counts_example.tsv")
#filter out zero counts to speed up algorithm
counts.dt<-counts.dt[rowSums(counts.dt[,.SD>0,.SD=c(2:dim(counts.dt)[2])])>0,]
probs4detection<-fread("tutorial_example_bayesianestimates.xls")
parameters<-fread("tutorial_example_parameters4randomize.xls")

#### Simulate replicates ####
cl <- makeCluster(no_cores, FORK=TRUE)
counts.error<-prepare_probabilities(counts.dt,
    probs4detection=probs4detection, parameters=parameters,
    HPC_genewise_permute_count=HPC_genewise_permute_count,
    HPC_permute_count=HPC_permute_count, HPC_randomizer=HPC_randomizer,
    total_sampling=2500)
counts.error.df<-data.frame(t(counts.error), row.names=counts.dt$GENE_ID)
counts.error.dt<-data.table(counts.error.df, keep.rownames=TRUE)
colnames(counts.error.dt)<-colnames(counts.dt)
write.table(counts.error.dt, file=paste("simulated_replicates/",
    paste(ITERATION,"sim_replicate_counts.txt",sep="_"),
    sep=""), quote =FALSE, row.names=FALSE)
stopCluster(cl)

#########
