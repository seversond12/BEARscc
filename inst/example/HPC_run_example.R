#!/usr/bin/env Rscript
#Author: David T Severson
#Scope: Example to run `prepare_probabilities()` on  a high perfromance
    #compute cluster with parallel performance

library("BEARscc")

#### Load data ####
ITERATION<-commandArgs(trailingOnly=TRUE)[1]
counts.df<-read.delim("tutorial_example_counts4clusterperturbation.xls")
#filter out zero counts to speed up algorithm
counts.df<-counts.df[rowSums(counts.df)>0,]

probs4detection<-fread("tutorial_example_bayesianestimates.xls")
parameters<-fread("tutorial_example_parameters4randomize.xls")

#### Simulate replicates ####
counts.error<-HPC_simulate_replicates(counts_matrix=counts.df,
    dropout_parameters=dropout_parameters,
    spikein_parameters=spikein_parameters)
write.table(counts.error, file=paste("simulated_replicates/",
    paste(ITERATION,"sim_replicate_counts.txt",sep="_"),
    sep=""), quote =FALSE, row.names=TRUE)
#########
