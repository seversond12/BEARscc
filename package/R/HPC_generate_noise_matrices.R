#!/usr/bin/env Rscript
#Author: David T Severson
#Scope: Noise injection matrices generation for use on high through-put clusters
library(data.table, logical.return=TRUE, quietly = TRUE)
library("parallel",logical.return=TRUE, quietly = TRUE)



#### function that permutes counts with technical noise, parameters are:
##NB_slope and NB_intercept are the result of relating the contribution of the negative binomial to the mean ERCC expression, where ##the contribution has the least residulas
##total_sampling is the number of points to sample in total split by n1 and n2 into negative binomial and poisson distributions
## Mean2SDcor_slope and Mean2SDcor_intercept are the result of fitting the standard deviation to the mean. HOwever, this squeezes
## the SD for ERCCs that are more dispersed than the other ERCCs as a funciton of mean are on the mean... 
##Mean2SDcor_residual describes the maximum residual observed in the SD to mean fit for ERCCs allowing higher dispersion
##Finally, the value x is simply the counts that lapply with use, and sd_inflate is a changeable parameter by the use that can be
## a real number from 0 to N, which allows the user to increase the spread by adding some multiple of the residual
##actually applies the estimated noise to the counts table
HPC_randomizer<-function(x,parameters,total_sampling){
  NB_intercept<-as.numeric(parameters$alpha2mu.intercept) 
  NB_slope<-as.numeric(parameters$alpha2mu.slope)
  Mean2SDcor_residual<-as.numeric(parameters$max.res)
  Mean2SDcor_slope<-as.numeric(parameters$mu2sigma.slope)
  Mean2SDcor_intercept<-as.numeric(parameters$mu2sigma.intercept)
  sd_inflate<-as.numeric(parameters$sd.inflate)
  p<-NB_intercept+NB_slope*log2(x+0.025)
  if (p<0) {p=0}
  if (p>1) {p=1}
  hx1 <- rnbinom(total_sampling*p,size = x^2/((2^(Mean2SDcor_slope*log2(x)+Mean2SDcor_intercept+Mean2SDcor_residual*sd_inflate))^2-x),mu = x)
  hx2 <- rpois(total_sampling*(1-p),x)
  sample(append(hx1, hx2), size=1, replace=TRUE)
}

##The count permuter that decides fate of observations of 0 and also calls the randomizer function.
HPC_permute_count<-function(x, probs4detection.k, probabilityA, parameters, gene, HPC_randomizer, total_sampling){
  force(parameters)
  if (x==0) {
    nx<-sample(as.numeric(rownames(probs4detection.k)),1,prob=probabilityA, replace=TRUE)
    if (nx>0) {
      nx<-HPC_randomizer(probs4detection.k[as.character(nx),]$Counts, parameters, total_sampling)
    }
  } 
  else {
    if (ceiling(x) < max(as.numeric(rownames(probs4detection.k)))){
      probabilityB<-as.numeric(probs4detection.k[as.character(ceiling(x)),]$Py0givenk)
      nx<-sample(c(0,x),1, prob=as.numeric(c(probabilityB, 1-probabilityB)))
      if (nx>0) {
        nx<-HPC_randomizer(nx, parameters, total_sampling)
      }
    } else {	
      nx<-HPC_randomizer(x, parameters, total_sampling)
    }
  }
  nx
}

##seperates out genes for analysis
HPC_genewise_permute_count<-function(x,probs4detection.k, probs4detection.genes, parameters,permute_count, randomizer,total_sampling=2500){
  print(x[1])
  print(sub("[-,(,)]",".",x[1]))
  probabilityA<-probs4detection.genes[gsub("[-,(,)]",".",x[1]),]
  force(total_sampling)
  apply(data.frame(as.numeric(x[-1])),1, `HPC_permute_count`, probs4detection.k,probabilityA, parameters,gene=gene_name,HPC_randomizer,total_sampling)
}

##prepares bayesian file for fast row indexing. 
prepare_probabilities<-function(x,probs4detection, parameters,HPC_genewise_permute_count=HPC_genewise_permute_count,HPC_permute_count=HPC_permute_count, HPC_randomizer=HPC_randomizer,total_sampling=2500){
  probs4detection.genes<-t(data.frame(probs4detection, row.names = "k")[,4:eval(dim(probs4detection)[2]-1)])
  rownames(probs4detection.genes)<-gsub("^X([0-9])","\\1", rownames(probs4detection.genes))
  probs4detection.k<-data.frame(probs4detection[,2:4, with=FALSE],row.names = "k")
  x[,parApply(cl,.SD, 1 ,`HPC_genewise_permute_count`,probs4detection.k=probs4detection.k,probs4detection.genes=probs4detection.genes,total_sampling=total_sampling, HPC_permute_count=HPC_permute_count, HPC_randomizer=HPC_randomizer, parameters=parameters)]
}

  


###Example template for HPC run of BEARscc noise counts matrix generation###
##Load data##
#no_cores<-4
#counts.dt<-fread("brain_control_example.tsv")
#counts.dt<-counts.dt[rowSums(counts.dt[,.SD>0,.SD=c(2:dim(counts.dt)[2])])>0,]
#probs4detection<-fread("example_bayesianestimates.xls")
#parameters<-fread("example_parameters4randomize.xls")
####

##error test##
#cl <- makeCluster(no_cores, FORK=TRUE)
#counts.error<-prepare_probabilities(counts.dt, probs4detection=probs4detection, parameters=parameters, HPC_genewise_permute_count=HPC_genewise_permute_count, HPC_permute_count=HPC_permute_count, HPC_randomizer=HPC_randomizer,total_sampling=2500)
#counts.error.df<-data.frame(t(counts.error), row.names=counts.dt$GENE_ID)
#counts.error.df<-round(counts.error.df)
#counts.error.dt<-data.table(counts.error.df, keep.rownames=TRUE)
#colnames(counts.error.dt)<-colnames(counts.dt)
#write.table(counts.error.dt, file=paste("noise_perturbed_count_tables/", paste(TASKID,"noise_perturbation","perturbed_counts.txt",sep="_"), sep=""), quote =FALSE, row.names=FALSE)
#stopCluster(cl)
#####



