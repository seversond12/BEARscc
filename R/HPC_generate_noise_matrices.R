#!/usr/bin/env Rscript
#Author: David T Severson
#Scope: Simulate replicates on high throughput cluster
library(data.table, logical.return=TRUE, quietly = TRUE)
library("parallel",logical.return=TRUE, quietly = TRUE)

##cluster number
cl<-2

##function that permutes counts with technical noise
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
    hx1 <- rnbinom(total_sampling*p,size = x^2/((2^(Mean2SDcor_slope*
        log2(x)+Mean2SDcor_intercept+Mean2SDcor_residual*
        sd_inflate))^2-x),mu = x)
    hx2 <- rpois(total_sampling*(1-p),x)
    sample(append(hx1, hx2), size=1, replace=TRUE)
}

##The count permuter permutes zeros and also calls randomizer()
HPC_permute_count<-function(x, probs4detection.k, probabilityA,
    parameters, gene, HPC_randomizer, total_sampling){
    force(parameters)
    if (x==0) {
        nx<-sample(as.numeric(rownames(probs4detection.k)), 1,
            prob=probabilityA, replace=TRUE)
        if (nx>0) {
            nx<-HPC_randomizer(probs4detection.k[as.character(nx),]$Counts,
                parameters, total_sampling)
        }
    }
    else {
        if (ceiling(x) < max(as.numeric(rownames(probs4detection.k)))){
            probabilityB<-as.numeric(
                probs4detection.k[as.character(ceiling(x)),]$Py0givenk)
            nx<-sample(c(0,x),1, prob=as.numeric(c(probabilityB,
                1-probabilityB)))
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
HPC_genewise_permute_count<-function(x,probs4detection.k,
    probs4detection.genes, parameters,permute_count, randomizer,
    total_sampling=2500){
    gene_name<-NULL
    print(x[1])
    print(sub("[-,(,)]",".",x[1]))
    probabilityA<-probs4detection.genes[gsub("[-,(,)]",".",x[1]),]
    force(total_sampling)
    apply(data.frame(as.numeric(x[-1])),1, `HPC_permute_count`,
        probs4detection.k,probabilityA, parameters,gene=gene_name,
        HPC_randomizer,total_sampling)
}

##prepares bayesian file for fast row indexing.
prepare_probabilities<-function(x,probs4detection, parameters,
    HPC_genewise_permute_count=HPC_genewise_permute_count,
    HPC_permute_count=HPC_permute_count, HPC_randomizer=HPC_randomizer,
    total_sampling=2500){
    probs4detection.genes<-t(data.frame(probs4detection,
        row.names = "k")[,4:eval(dim(probs4detection)[2]-1)])
    rownames(probs4detection.genes)<-gsub("^X([0-9])","\\1",
        rownames(probs4detection.genes))
    probs4detection.k<-data.frame(probs4detection[,2:4,
        with=FALSE],row.names = "k")
    x[,parApply(cl,.SD, 1 ,`HPC_genewise_permute_count`,
        probs4detection.k=probs4detection.k,
        probs4detection.genes=probs4detection.genes,
        total_sampling=total_sampling, HPC_permute_count=HPC_permute_count,
        HPC_randomizer=HPC_randomizer, parameters=parameters)]
}
