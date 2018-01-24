#!/usr/bin/env Rscript
#Author: David T Severson
#Scope: Simulate replicates on high throughput cluster
library("data.table", quietly = TRUE)

HPC_simulate_replicates<-function(counts_matrix, dropout_parameters,
spikein_parameters, max_cumprob=0.9999){
    round_counts<-counts<-transcripts<-dropouts_given_transcripts<-dif<-NULL
    max_cumprob<-1-(1-max_cumprob)/2
    dropout_recovery.genes<-t(data.frame(dropout_parameters,
        row.names = "transcripts")[,3:eval(dim(
        dropout_parameters)[2]-1)])
    rownames(dropout_recovery.genes)<-gsub("^X([0-9])","\\1",
        rownames(dropout_recovery.genes))
    dropout_recovery<-data.table(dropout_parameters[,1:3, with=FALSE])[,
        round_counts:=round(counts),]
    dropout_injection<-data.table(dropout_parameters[,1:3, with=FALSE])[,
        round_counts:=round(counts),][,`:=`(transcripts=transcripts,
        dropouts_given_transcripts=dropouts_given_transcripts,
        dif=abs(counts-round_counts), min=min(abs(counts-round_counts))),
        by=round_counts][dif==min,][,`:=`(transcripts=transcripts, dif=NULL,
        min=NULL, dropouts_given_transcripts=dropouts_given_transcripts,
        counts=round_counts, round_counts=NULL)]
    count_max<-max(dropout_injection$counts)
    all_counts<-seq(from=1,to=max(dropout_injection$counts), by=1)
    if (length(all_counts)>(length(dropout_injection$counts)-1)){
        dropout_injection<-do.call(rbind, lapply(all_counts,
            `fill_out_count_probability_table`, dropout_injection))
        dropout_injection<-data.frame(dropout_injection, row.names="counts")
    } else{
        dropout_injection<-data.frame(dropout_injection, row.names = "counts")
    }
    noisy_counts<-t(counts_matrix)
    noisy_counts.zeros<-noisy_counts[,colSums(noisy_counts==0)>0]
    noisy_counts.nozeros<-noisy_counts[,colSums(noisy_counts==0)==0]
    indropout_range<-which(noisy_counts<=count_max & noisy_counts>0)
    dropouts<-sapply(noisy_counts[indropout_range],
        `permute_count_in_dropout_range`, dropout_injection)
    t_counts<-rbind(colnames(noisy_counts.zeros), noisy_counts.zeros)
    noisy_counts.zeros<-data.table(t_counts)[,sapply(.SD,
        `genewise_dropouts`, dropout_recovery=dropout_recovery,
        dropout_recovery.genes=dropout_recovery.genes)]
    noisy_counts<-cbind(noisy_counts.zeros, noisy_counts.nozeros)[
        rownames(noisy_counts), colnames(noisy_counts)]
    noisy_counts[indropout_range]<-dropouts
    noisy_counts[noisy_counts>0]<-sapply(noisy_counts[noisy_counts>0],
        `randomizer`, parameters=spikein_parameters, max_cumprob=max_cumprob)
    t(noisy_counts)
}
