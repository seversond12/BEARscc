# R script
# David Severson
# Scope: To simulate technical replicates.

################Underlying functions############################

##actually applies the estimated noise to the counts table
randomizer<-function(count, parameters, max_cumprob){
    qnbinom<-qpois<-dnbinom<-dpois<-NULL
    count<-ceiling(count)
    NB_intercept<-as.numeric(parameters$alpha2mu.intercept)
    NB_slope<-as.numeric(parameters$alpha2mu.slope)
    Mean2SDcor_residual<-as.numeric(parameters$max.res)
    Mean2SDcor_slope<-as.numeric(parameters$mu2sigma.slope)
    Mean2SDcor_intercept<-as.numeric(parameters$mu2sigma.intercept)
    sd_inflate<-as.numeric(parameters$sd.inflate)
    p<-NB_intercept+NB_slope*log2(count)
    if (p<0) {p=0}
    if (p>1) {p=1}
    top_probable_count<-max(qnbinom(max_cumprob, size = count^2/
        ((2^(Mean2SDcor_slope*log2(count)+Mean2SDcor_intercept+
        Mean2SDcor_residual*sd_inflate))^2-count),mu = count),
        qpois(max_cumprob, count))
    min_probable_count<-min(qnbinom(max_cumprob, size = count^2/
        ((2^(Mean2SDcor_slope*log2(count)+Mean2SDcor_intercept+
        Mean2SDcor_residual*sd_inflate))^2-count), mu = count,
        lower.tail = FALSE), qpois(max_cumprob, count, lower.tail = FALSE))
    sequence<-seq(from=min_probable_count, to=top_probable_count, by=1)
    sample(sequence,size=1,prob=p*dnbinom(sequence, size = count^2/((2^(
        Mean2SDcor_slope*log2(count)+Mean2SDcor_intercept+
        Mean2SDcor_residual*sd_inflate))^2-count), mu = count)+(1-p)*dpois(
        sequence,count))
}

permute_count_in_dropout_range<-function(count, dropout_injection){
    count<-round(count)
    dropout_prob<-dropout_injection[count,]$dropouts_given_transcripts
    sample(c(0,count),1, prob=as.numeric(c(dropout_prob, 1-dropout_prob)))
}

##seperates out genes for analysis
genewise_dropouts<-function(onegene_counts, dropout_recovery,
dropout_recovery.genes){
    round_counts<-NULL
    genewise_dropout_recovery_prob<-dropout_recovery.genes[gsub("[-,(,)]",
        ".",onegene_counts[1]),]
    onegene_counts<-as.numeric(onegene_counts[-1])
    sampling<-sample(
      dropout_recovery$transcripts, sum(onegene_counts==0),
      prob=genewise_dropout_recovery_prob, replace=TRUE)
    onegene_counts[onegene_counts==0]<-dropout_recovery[match(sample(
        dropout_recovery$transcripts, sum(onegene_counts==0),
        prob=genewise_dropout_recovery_prob, replace=TRUE),
        dropout_recovery$transcripts), round_counts]
    onegene_counts
}

fill_out_count_probability_table<-function(count, dropout_injection){
    counts<-transcripts<-NULL
    dropout_injection[counts>=count,][,min:=min(transcripts),][
      transcripts==min,][,counts:=count][,min:=NULL]
}

execute_sim_replicates<-function(sim_replicate, SCEList, max_cumprob){
    metadata<-round_counts<-counts<-transcripts<-NULL
    dropouts_given_transcripts<-dif<-assay<-NULL
    print(paste("Creating a simulated replicated counts matrix: ",
        sim_replicate,".", sep=""))
    noise_parameters<-metadata(SCEList)
    dropout_recovery.genes<-t(data.frame(noise_parameters$dropout_parameters,
        row.names = "transcripts")[,3:eval(dim(
        noise_parameters$dropout_parameters)[2]-1)])
    rownames(dropout_recovery.genes)<-gsub("^X([0-9])","\\1",
        rownames(dropout_recovery.genes))
    dropout_recovery<-data.table(noise_parameters$dropout_parameters[,1:3])[,
        round_counts:=round(counts),]
    dropout_injection<-data.table(noise_parameters$dropout_parameters[,1:3])[,
        round_counts:=round(counts),][,`:=`(transcripts=transcripts,
        dropouts_given_transcripts=dropouts_given_transcripts,
        dif=abs(counts-round_counts), min=min(abs(counts-round_counts))),
        by=round_counts][dif==min,][,`:=`(transcripts=transcripts, dif=NULL,
        min=NULL, dropouts_given_transcripts=dropouts_given_transcripts,
        counts=round_counts, round_counts=NULL)]
    count_max<-max(dropout_injection$counts)
    all_counts<-seq(from=1,to=max(dropout_injection$counts), by=1)
    if (length(all_counts)>(length(dropout_injection$counts)-1)){
        dropout_injection<-data.frame(do.call(rbind, lapply(all_counts,
           `fill_out_count_probability_table`, dropout_injection)),
           row.names="counts")
    } else{
        dropout_injection<-data.frame(dropout_injection, row.names = "counts")
    }
    noisy_counts<-t(data.frame(assay(SCEList, "observed_expression")))
    noisy_counts.zeros<-noisy_counts[,colSums(noisy_counts==0)>0]
    noisy_counts.nozeros<-noisy_counts[,colSums(noisy_counts==0)==0]
    indropout_range<-which(noisy_counts.zeros<=count_max & noisy_counts.zeros>0)
    dropouts<-sapply(noisy_counts.zeros[indropout_range],
        `permute_count_in_dropout_range`, dropout_injection)
    t_counts<-rbind(colnames(noisy_counts.zeros), noisy_counts.zeros)
    noisy_counts.zeros<-data.table(t_counts)[,sapply(.SD,`genewise_dropouts`,
        dropout_recovery=dropout_recovery,
        dropout_recovery.genes=dropout_recovery.genes)]
    noisy_counts<-cbind(noisy_counts.zeros, noisy_counts.nozeros)[
        rownames(noisy_counts), colnames(noisy_counts)]
    noisy_counts.zeros[indropout_range]<-dropouts
    noisy_counts[noisy_counts>0]<-sapply(noisy_counts[noisy_counts>0],
        `randomizer`, parameters=noise_parameters$spikein_parameters,
         max_cumprob=max_cumprob)
    t(noisy_counts)
}

############# USER FUNCTIONS ###################################

##executes noise perturbation with estimated parameters
simulate_replicates<-function(SCEList, max_cumprob=0.9999, n=3){
    metadata<-`metadata<-`<-NULL
    max_cumprob<-1-(1-max_cumprob)/2
    sim_rep_vector<-data.table(t(as.character(seq(from=1, to=n, by=1))))
    colnames(sim_rep_vector)<-gsub("^","Iteration_", sim_rep_vector)
    metadata(SCEList)$simulated_replicates<-lapply(sim_rep_vector,
        `execute_sim_replicates`, SCEList=SCEList, max_cumprob=max_cumprob)
    SCEList
}


