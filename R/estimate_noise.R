#! R script
# David T Severson
# Scope: Function in BEARscc used to estimate parameters
    # and dropouts from ERCC spike-ins.


################### INTERNAL FUNCTIONS #######################################


#annotate spike-in counts with actual transcript counts
prepare_data<-function(spike_counts.df,spike_conc.df){
    spike_counts.df$transcripts<-spike_conc.df[rownames(spike_counts.df),1]
    spike_counts.df
}

#melt and format spike-in counts
melt_spikeins<-function(prepared_spikeins){
    value<-V1<-NULL
    melted_spikeins.dt<-melt(data.table(prepared_spikeins,
        keep.rownames = TRUE), id.vars = c("rn", "transcripts"))
    detected_spikeins.dt<-melted_spikeins.dt[,mean(value)>0, by="rn"]
    setkey(detected_spikeins.dt, "rn")
    setkey(melted_spikeins.dt,"rn")
    melted_spikeins.dt<-detected_spikeins.dt[melted_spikeins.dt,
        ][V1==TRUE,][,-2,with=FALSE]
    colnames(melted_spikeins.dt)<-c("spike-in", "transcripts" , "sample",
        "counts")
    melted_spikeins.dt
}

#plot mean to std dev correlation
plot_mu2sigma<-function(spikein_stats.dt, mu2sigma.fit, file, sd_inflate){
    V1<-NULL
    g<-ggplot(data=spikein_stats.dt,aes(x=log2(mean), y=log2(V1)))+
        geom_point()+xlab("\nMean spike-in expression / log2")+
        geom_abline(slope = coef(mu2sigma.fit)[2],
        intercept = coef(mu2sigma.fit)[1])+
        ylab("Standard deviation of spike-in expression / log2\n")+
    ggtitle("Relationship between standard deviation
        and mean of spike in measurements\n")+
    geom_abline(slope=coef(mu2sigma.fit)[2],color="red",
        intercept=coef(mu2sigma.fit)[1]+
        sd_inflate*max(mu2sigma.fit$residuals))
    ggsave(paste(file,"mu2sigma.pdf", sep="_"), plot=g, device="pdf",
        width=8, height=6, units = "in")
}

#estimate mean to std dev correlation by linear regression
estimate_mu2sigma<-function(melted_spikeins.dt, plot, sd_inflate, file=file){
    V1<-value<-counts<-NULL
    spikein_stats.dt<-melted_spikeins.dt[,sd(counts), by=c("spike-in",
        "transcripts")]
    spikein_stats.dt$mean<-melted_spikeins.dt[,mean(counts),
        by=c("spike-in", "transcripts")]$V1
    mu2sigma.fit<-lm(log2(spikein_stats.dt$V1)~log2(spikein_stats.dt$mean))
    mu2sigma<-data.frame(mu2sigma.slope=coef(mu2sigma.fit)[2],
        mu2sigma.intercept=coef(mu2sigma.fit)[1],
        max.res=max(mu2sigma.fit$residuals))
    if (plot==TRUE){
        plot_mu2sigma(spikein_stats.dt, mu2sigma.fit, file=file,
            sd_inflate=sd_inflate)
    }
    mu2sigma
}

#plot actual spike-in transcript number vs observed spike-in counts
plot_obs2actual<-function(melted_spikeins.dt, file){
    transcripts<-counts<-NULL
    g<-ggplot(data=melted_spikeins.dt[transcripts>0,],
        aes(x=factor(round(log10(transcripts+0.25), digits=3)),
        y=log2(counts+0.25)))+geom_violin(scale = "width")+
    xlab("\nSpike-in molecules per cell / log10")+
    ylab("Measured counts / log2\n")+
    ggtitle("The count measurement distribution as a
        function of actual transcript number.\n")
    ggsave(paste(file,"_countsdistribution_vs_moleculespercell.pdf", sep="_"),
        plot=g, device="pdf", width=8, height=6, units = "in")
}

#subfunction to compute alpha_i for each spike-in
iterate_spikeins<-function(spikein_name, alpha, melted_spikeins.dt,
mean_spikeins.dt, max_cumprob, bins, mu2sigma, sd_inflate){
    `spike-in`<-NULL
    melted_1spikein.dt<-melted_spikeins.dt[`spike-in`==spikein_name]
    onespikein_samples<-dim(melted_1spikein.dt)[1]
    mean_1spikein<-mean_spikeins.dt[`spike-in`==spikein_name]$mean
    top_probable_count<-max(qnbinom(max_cumprob,
        size=ceiling(mean_1spikein)^2/((2^(mu2sigma[1,1]*
        log2(ceiling(mean_1spikein))+mu2sigma[1,2]+
        sd_inflate*mu2sigma[1,3]))^2-ceiling(mean_1spikein)),
        mu=ceiling(mean_1spikein)), qpois(max_cumprob,
        ceiling(mean_1spikein)), melted_1spikein.dt$counts+1,
        na.rm = TRUE)
    event_space<-seq(from=0, to=top_probable_count, by=1)
    neg_binom.density <- dnbinom(event_space,
        size=ceiling(mean_1spikein)^2/((2^(mu2sigma[1,1]*
        log2(ceiling(mean_1spikein))+mu2sigma[1,2]+
        sd_inflate*mu2sigma[1,3]))^2-round(mean_1spikein)),
        mu=ceiling(mean_1spikein))
    poisson.density <- dpois(event_space,ceiling(mean_1spikein))
    mixed_model.density<-alpha*neg_binom.density+(1-alpha)*poisson.density
    event_space2<-seq(from=0, to=top_probable_count,
        by=ceiling(top_probable_count/bins))
    observed.hist<-hist(melted_1spikein.dt$counts, breaks=event_space,
        include.lowest = TRUE, plot=FALSE)
    observed.density<-observed.hist$density
    observed.density<-diff(c(0,cumsum(observed.density)[event_space2]))
    mixed_model.density<-diff(c(0, cumsum(mixed_model.density)[event_space2]))
    error<-sum(abs(mixed_model.density-observed.density))
    data.frame(`spike-in`=as.character(spikein_name),mean=mean_1spikein,
    alpha=alpha, alpha_inv=1-alpha, error=error)
}

#subfunction to compute alpha_i for each spike-in
iterate_alphas<-function(alpha, spikein_names, melted_spikeins.dt,
mean_spikeins.dt, max_cumprob, bins, mu2sigma, sd_inflate){
    if (alpha %in% seq(from=0, to=1, by=0.25)){
        print(paste("Estimating error for spike-ins with alpha = ", alpha,
            sep=""))
    }
    alpha_table.df<-do.call(rbind, lapply(spikein_names, `iterate_spikeins`,
        alpha, melted_spikeins.dt, mean_spikeins.dt, max_cumprob, bins,
        mu2sigma, sd_inflate))
    alpha_table.df
}

#plots the alpha paramter as function of mean spike-in expression
plot_alpha2mu<-function(final_alpha.dt, alpha2mean, file){
    g<-ggplot(data=final_alpha.dt, aes(x=log2(mean),y=alpha))+
        geom_point(fill="red",color="black", size=2.2,alpha=0.6,  pch=21)+
        xlab("\nObserved mean spike-in expression, log2(expression)")+
        ylab("Neg. binomial contribution parameter, alpha\n")+ggtitle("The Neg.
        binomial contribution is a function of gene expression.\n")+
        geom_abline(slope=coef(alpha2mean)[2], intercept = coef(alpha2mean)[1])
    ggsave(paste(file,"alpha2mu_correlation.pdf",sep="_"),
        plot=g, device="pdf", width=8, height=6, units = "in")
}

#compute alpha parameter empirically for mixed-model
compute_alpha<-function(melted_spikeins.dt,estimate_mu2sigma, bins,
tie_function, max_cumprob, sd_inflate, plot, file,
alpha_resolution=alpha_resolution){
    transcripts<-value<-rn<-Mean<-r.value<-p.nbinom<-NULL
    counts<-error<-`spike-in`<-min_error<-NULL
    mu2sigma<-estimate_mu2sigma(melted_spikeins.dt, plot, sd_inflate, file)
    if (plot==TRUE){
        plot_obs2actual(melted_spikeins.dt, file=file)
    }
    mean_spikeins.dt<-melted_spikeins.dt[,mean(counts),
        by=c("spike-in", "transcripts")]
    setkey(mean_spikeins.dt, 'spike-in')
    colnames(mean_spikeins.dt)<-c('spike-in', "transcripts", "mean")
    spikein_names<-as.character(mean_spikeins.dt$'spike-in')
    setkey(mean_spikeins.dt, 'spike-in')
    alpha_vector<-seq(from=0, to=1, by=alpha_resolution)
    alpha_table.dt<-data.table(do.call(rbind, lapply(alpha_vector,
        `iterate_alphas`, spikein_names, melted_spikeins.dt,
        mean_spikeins.dt, max_cumprob, bins, mu2sigma, sd_inflate)))
    colnames(alpha_table.dt)<-c("spike-in", "mean", "alpha", "alpha_inv",
        "error")
    min_error.dt<-alpha_table.dt[,min(error),by=`spike-in`]
    colnames(min_error.dt)<-c("spike-in", "min_error")
    setkey(min_error.dt, `spike-in`)
    setkey(alpha_table.dt, `spike-in`)
    alpha_table.dt<-alpha_table.dt[min_error.dt]
    if (tie_function=="minimum") {
        final_alpha.dt<-alpha_table.dt[error==min_error, min(alpha),
        by=c("spike-in", "mean")]
    }
    else {
        final_alpha.dt<-alpha_table.dt[error==min_error, min(alpha),
        by=c("spike-in", "mean")]
    }
    colnames(final_alpha.dt)<-c("spike-in", "mean", "alpha")
    alpha2mean<-lm(final_alpha.dt$alpha~log2(final_alpha.dt$mean))
    if (plot==TRUE){
        plot_alpha2mu(final_alpha.dt, alpha2mean, file=file)
    }
    alpha2mu<-data.frame(alpha2mu.slope=coef(alpha2mean)[2],
        alpha2mu.intercept=coef(alpha2mean)[1], sd.inflate=sd_inflate)
    model_parameters<-data.frame(cbind(mu2sigma,alpha2mu),
        row.names=c("parameter.value"))
    model_parameters
}

#spike-in wise model computation
subcompute_sample_models<-function(spikein_name, melted_spikeins.dt,
model_parameters){
    `spike-in`<-NULL
    one_spikein.dt<-melted_spikeins.dt[`spike-in`==spikein_name]
    transcript_conc<-melted_spikeins.dt[`spike-in`==spikein_name]$transcripts
    spikenamerepeats<-melted_spikeins.dt[`spike-in`==spikein_name]$`spike-in`
    sample_number<-length(one_spikein.dt$counts)
    mean_1spikein<-mean(one_spikein.dt$counts)
    alpha<-model_parameters$alpha2mu.slope*log2(mean_1spikein)+
        model_parameters$alpha2mu.intercept
    if (alpha>1) {alpha<-1}
    if (alpha<0) {alpha<-0}
    neg_binom.sample <- rnbinom(sample_number, size = ceiling(mean_1spikein
        )^2/((2^(model_parameters$mu2sigma.slope*log2(ceiling(mean_1spikein))+
        model_parameters$mu2sigma.intercept+model_parameters$max.res*
        model_parameters$sd.inflate))^2-ceiling(mean_1spikein)),
        mu = ceiling(mean_1spikein))
    poisson.sample <- rpois(sample_number, mean_1spikein)
    neg_binom_mixed.sample <- rnbinom(ceiling(sample_number*alpha),
        size = ceiling(mean_1spikein)^2/((2^(model_parameters$mu2sigma.slope*
        log2(ceiling(mean_1spikein))+model_parameters$mu2sigma.intercept+
        model_parameters$max.res*model_parameters$sd.inflate
        ))^2-ceiling(mean_1spikein)), mu = ceiling(mean_1spikein))
    poisson_mixed.sample <- rpois(ceiling(sample_number*(1-alpha)),
        mean_1spikein)
    mixed_model.sample<-sample(append(neg_binom_mixed.sample,
        poisson_mixed.sample), sample_number, replace = FALSE)
    samp.mixedmodel<-data.frame(`spike-in`=spikenamerepeats,
        transcripts=as.numeric(transcript_conc), sample=rep("Mixed",
        sample_number), counts=mixed_model.sample)
    samp.negbinomial<-data.frame(`spike-in`=spikenamerepeats,
        transcripts=transcript_conc, sample=rep("NBinom", sample_number),
        counts=neg_binom.sample)
    samp.poisson<-data.frame(`spike-in`=spikenamerepeats,
        transcripts=transcript_conc, sample=rep("Poisson", sample_number),
        counts=poisson.sample)
    model_sampling.df<-rbind(samp.mixedmodel, samp.negbinomial,
        samp.poisson, data.frame(one_spikein.dt))
    model_sampling.df
}

#computes vanilla, and mixed models for plotting
sample_models<-function(melted_spikeins.dt, model_parameters){
    rn<-NULL
    spikein_names<-unique(as.character(melted_spikeins.dt$`spike-in`))
    model_sampling<-do.call(rbind,lapply(spikein_names,
        `subcompute_sample_models`, melted_spikeins.dt, model_parameters))
    data.table(model_sampling)
}

#makes plots of models more readable
cleanup_model_names<-function(model_sampling){
    model_sampling$distribution<-gsub(TRUE,"Observed",
        model_sampling$sample!="Mixed" & model_sampling$sample!="NBinom" &
        model_sampling$sample!="Poisson")
    model_sampling[model_sampling$sample=="Mixed",
        ]$distribution<-"Optimized"
    model_sampling[model_sampling$sample=="Poisson",
        ]$distribution<-"Poisson"
    model_sampling[model_sampling$sample=="NBinom",
        ]$distribution<-"Neg. Binomial"
    data.table(model_sampling)
}

subplot_spikein_fits<-function(iteration,model_sampling.dt, model_view, file,
bins, spikein_names){
    spike.in<-distribution<-counts<-NULL
    six_spikeins<-spikein_names[seq(from=iteration, to=length(spikein_names),
        by=round(length(spikein_names)/6))]
    g<-ggplot(data=model_sampling.dt[spike.in %in% six_spikeins,
        ][distribution %in% model_view], aes(x=counts))+
        geom_histogram(aes(fill=factor(distribution, levels=c("Poisson",
        "Neg. Binomial", "Optimized", "Observed"))), alpha=0.30,
        position="identity", bins=bins)+facet_wrap(~spike.in, ncol=3,nrow=4,
        scales="free")+ylab("Observation density\n")+labs(fill="Model")+
        xlab("\nSpike-in expression / normalized counts")+
        ggtitle("Approximating technical noise using a refined
        mixed distribution.\n")+scale_color_discrete(guide=FALSE)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10))
    savefile<-paste("spikein_fits/", paste(file,"spikein_fits_subset",
        iteration, "plot.pdf",sep="_"),sep="")
    ggsave(savefile, plot=g, device="pdf", width=8, height=6, units = "in")
}

#plot mixed-model fits for each spike-in molecule
plot_spikein_fits<-function(model_sampling.dt, model_view, file, bins){
    transcripts<-NULL
    dir.create("spikein_fits")
    print("Writing spike-in plots; this takes a few seconds.")
    spikein_names<-as.character(unique(model_sampling.dt[
        order(transcripts)]$spike.in))
    iteration<-seq(from=1, to=round(length(spikein_names)/6), by=1)
    lapply(iteration, `subplot_spikein_fits`,
        model_sampling.dt=model_sampling.dt, model_view=model_view,
        file=file, bins=bins, spikein_names=spikein_names)
}

#computes the linear fit parameters to model false zeroes
estimate_undetected2molpercell<-function(undetected_spikeins.dt,
plot, file){
    V1<-transcripts<-undetected<-NULL
    undetected_spikeins.dt<-undetected_spikeins.dt[
        undetected!=0][transcripts>=0.5]
    undetected2mpc.fit<-lm(undetected~log2(transcripts),
        undetected_spikeins.dt)
    ##Read correlation paramaters into data.frame
    undetected2molpercell<-data.frame(
        undetected2mpc.slope=coef(undetected2mpc.fit)[2],
        undetected2mpc.intercept=coef(undetected2mpc.fit)[1])
    if (plot==TRUE){
        g<-ggplot(data=undetected_spikeins.dt,
            aes(x=log2(transcripts), y=undetected))+geom_point()+
            xlab("\nSpike-in actual counts, log2(molecules per cell)")+
            ylab("Fraction of cells where observed spike-in count is zero")+
            ggtitle("Conditional probability of observing a zero given
            transcript concentration\n")+geom_abline(
            slope=coef(undetected2mpc.fit)[2],color="red",
            intercept=coef(undetected2mpc.fit)[1])
        ggsave(paste(file,"undetected2molpercell.pdf", sep="_"), plot=g,
            device="pdf", width=11, height=8.5, units = "in")
    }
    undetected2molpercell
}

#computes relationship between counts and molecules per cell
counts2mpc<-function(melted_spikeins.dt, plot, file){
    transcripts<-value<-V1<-counts<-NULL
    mean_spikeins.dt<-melted_spikeins.dt[,mean(counts),by=transcripts][
        transcripts>0]
    counts2mpc.fit<-lm(log2(mean_spikeins.dt$V1)~log2(
        mean_spikeins.dt$transcripts))
    if (plot==TRUE){
        g<-ggplot(data=mean_spikeins.dt, aes(x=log2(transcripts),y=log2(V1)))+
            geom_point()+xlab("\nSpike-in molecules per cell,
            log2(molecules)")+ylab("Observed counts, log2(counts) \n")+
            ggtitle("Correlation of measured counts as afunction of actual
            transcript number.\n")+geom_abline(slope=coef(counts2mpc.fit)[2],
            intercept=coef(counts2mpc.fit)[1], color="red")
        ggsave(paste(file,"cor_counts2moleculespercell.pdf", sep="_"),
            plot=g, device="pdf", width=8, height=6, units = "in")
    }
    counts2mpc.fit
}

apply_bayes<-function(dropout_model.dt, noise_models){
    dropouts_given_transcripts<-as.vector(
        dropout_model.dt$dropouts_given_transcripts[-1])
    genewise_prob_transcript_is_k<-as.vector(
        noise_models$genewise$prob_transcript_is_k)
    names(genewise_prob_transcript_is_k)<-as.character(
        rownames(noise_models$genewise_dropouts))
    bayes_numerator<-dropouts_given_transcripts%*%t(
        genewise_prob_transcript_is_k)
    noise_models$dropout_parameters<-data.table(data.frame(sweep(rbind(t(
        noise_models$genewise_dropouts$prob_truezero), bayes_numerator), 2,
        noise_models$genewise_dropouts$prob_geneobservedaszero, "/")),
        keep.rownames=TRUE)
    colnames(noise_models$dropout_parameters)<-c("transcripts",
        colnames(noise_models$dropout_parameters)[-1])
    noise_models$dropout_parameters$transcripts<-
        as.numeric(noise_models$dropout_parameters$transcripts)-1
    noise_models
}

#computes genewise probabilities that observed values are zero
compute_genewise_dropouts<-function(dropout_model.dt, endogenous_counts.df,
kmax){
    transcripts<-NULL
    noise_models<-list()
    noise_models$genewise_dropouts<-data.frame(rowSums(endogenous_counts.df==0
        )/dim(endogenous_counts.df)[2])
    colnames(noise_models$genewise_dropouts)<-"prob_geneobservedaszero"
    noise_models$genewise_dropouts$prob_truezero<-(
        noise_models$genewise_dropouts$prob_geneobservedaszero-
        sum(dropout_model.dt$dropouts_given_transcripts[-1])/kmax)/
        (1-sum(dropout_model.dt$dropouts_given_transcripts[-1])/kmax)
    noise_models$genewise_dropouts[
        noise_models$genewise_dropouts$prob_truezero<0,]$prob_truezero<-0
    noise_models$genewise_dropouts$prob_transcript_is_k<-(1-
        noise_models$genewise_dropouts$prob_truezero)*1/kmax
    noise_models<-apply_bayes(dropout_model.dt=dropout_model.dt,
        noise_models=noise_models)
    setkey(noise_models$dropout_parameters, transcripts)
    setkey(dropout_model.dt, transcripts)
    noise_models$dropout_parameters<-dropout_model.dt[
        noise_models$dropout_parameters,]
    setkey(noise_models$dropout_parameters, transcripts)
    noise_models$dropout_parameters<-data.frame(
        noise_models$dropout_parameters)
    noise_models
}

build_dropoutmodel<-function(kmax, undetected2mpc, counts2mpc.fit,
dropout_inflate){
    transcript_values<-seq(from=0, to=kmax, by=1)
    dropout_model.dt<-data.table(data.frame(
        transcripts=transcript_values))
    dropout_model.dt$dropouts_given_transcripts<-(1/as.numeric(
        dropout_inflate))*as.numeric(undetected2mpc[1])*log2(
        dropout_model.dt$transcripts)+as.numeric(undetected2mpc[2])
    dropout_model.dt[dropout_model.dt$dropouts_given_transcripts>1 |
        dropout_model.dt$dropouts_given_transcripts==Inf,
        ]$dropouts_given_transcripts<-1
    dropout_model.dt$counts<-2^(coef(counts2mpc.fit)[2]*
        log2(dropout_model.dt$transcripts)+coef(counts2mpc.fit)[1])
    dropout_model.dt
}


create_null_dropout_model<-function(counts2mpc.fit, endogenous_counts.df){
    noise_models<-list()
    transcripts<-seq(from=0, to=10, by=1)
    dropouts_given_transcripts<-c(1,rep(0,10))
    counts<-c(0, 2^(coef(counts2mpc.fit)[2]*
        log2(transcripts[-1])+coef(counts2mpc.fit)[1]))
    gene_nodropouts<-data.frame(dropouts_given_transcripts%*%t(rep(1,
        length(rownames(endogenous_counts.df)))), row.names = transcripts)
    colnames(gene_nodropouts)<-rownames(endogenous_counts.df)
    dropout_null<-data.frame(transcripts=transcripts,
        dropouts_given_transcripts=dropouts_given_transcripts,
        counts=counts)
    rownames(dropout_null)<-transcripts
    noise_models$dropout_parameters<-cbind(dropout_null, gene_nodropouts)
    noise_models
}

#Uses Bayes' theorem to build drop-out noise models
estimate_missingdata<-function(melted_spikeins.dt, endogenous_counts.df,
counts2mpc.fit, plot, file, dropout_inflate=dropout_inflate, sample_number){
    transcripts<-rn<-counts<-NULL
    force(sample_number)
    undetected_spikeins.dt<-melted_spikeins.dt[,sum(counts==0)/sample_number,
        by=c("spike-in", "transcripts")][transcripts>0,]
    colnames(undetected_spikeins.dt)<-c("spike-in", "transcripts",
        "undetected")
    undetected2mpc<-estimate_undetected2molpercell(undetected_spikeins.dt,
        plot=plot, file=file)
    kmax<-ceiling(2^(-as.numeric(undetected2mpc[2])/
        ((1/as.numeric(dropout_inflate))*as.numeric(undetected2mpc[1]))))-1
    if (sum(undetected_spikeins.dt$undetected==0)==0){
        print(paste("Warning: there are no spike-ins that were detected in",
            "every sample. As a result the actual transcript count",
            "threshold, k, at which drop-outs are not present will be",
            "extrapolated rather than interpolated. The extrapolated value ",
            "for k is, ", kmax, ".", sep=""))
    }
    if (kmax>0){
        print(paste("There are adequate spike-in drop-outs to build the",
            "drop-out model. Estimating the drop-out model now."), sep="")
        noise_models<-compute_genewise_dropouts(build_dropoutmodel(kmax=kmax,
            undetected2mpc=undetected2mpc, counts2mpc.fit=counts2mpc.fit,
            dropout_inflate=dropout_inflate), kmax=kmax,
            endogenous_counts.df=endogenous_counts.df)
    } else {
        print(paste("Warning: there are inadequate spike-in drop-outs to",
            "build the drop-out model. All observed zeros will be treated",
            "as true zeros during simulation of technical replicates."),
            sep="")
        noise_models<-create_null_dropout_model(counts2mpc.fit=counts2mpc.fit,
            endogenous_counts.df=endogenous_counts.df)
    }
    noise_models
}

write_noise_model<-function(SCEList, file){
    metadata<-assay<-NULL
    write.table(metadata(SCEList)$spikein_parameters, file=paste(file,
        "parameters4randomize.xls", sep="_"), quote=FALSE,
        sep="\t", row.names = FALSE)
    write.table(metadata(SCEList)$dropout_parameters,
        file=paste(file, "bayesianestimates.xls", sep="_"),
        quote=FALSE, sep="\t", row.names=FALSE )
    write.table(assay(SCEList, "observed_expression"), file=paste(file,
        "counts4clusterperturbation.xls", sep="_"), quote=FALSE,
        sep="\t", row.names=TRUE)
    print(paste("Parameters have been saved into the current directory",
        "in bayesianestimates.xls and parameters4randomize.xls."), sep="")
}


############################ USER COMMANDS ###################################
estimate_noiseparameters<-function(SCEList, plot=FALSE, sd_inflate=0,
max_cumprob=0.9999, bins=10, write.noise.model=TRUE, file="noise_estimation",
dropout_inflate=1, model_view=c("Observed", "Optimized"),
alpha_resolution=0.005, tie_function="maximum"){
    transcripts<-value<-V1<-rn<-distribution<-metadata<-NULL
    `metadata<-`<-assay<-NULL
    spike_counts.df<-data.frame(assay(SCEList, "observed_expression")[isSpike(
        SCEList, "ERCC_spikes"),])
    spike_conc.df<-metadata(SCEList)$spikeConcentrations
    endogenous_counts.df<-data.frame(assay(SCEList, "observed_expression"))
    melted_spikeins.dt<-melt_spikeins(prepare_data(spike_counts.df,
        spike_conc.df))
    print("Fitting parameter alpha to establish spike-in derived noise model.")
    model_parameters<-compute_alpha(melted_spikeins.dt=melted_spikeins.dt,
        estimate_mu2sigma, bins=bins, max_cumprob=max_cumprob, plot=plot,
        sd_inflate=sd_inflate, file=file, tie_function = tie_function,
        alpha_resolution=alpha_resolution)
    model_sampling.dt<-cleanup_model_names(sample_models(melted_spikeins.dt,
        model_parameters))
    if (plot==TRUE){
        plot_spikein_fits(model_sampling.dt, model_view=model_view,
            file=file, bins=bins)
    }
    counts2mpc.fit<-counts2mpc(melted_spikeins.dt, plot=plot, file=file)
    sample_number<-dim(spike_counts.df)[2]
    noise_models<-estimate_missingdata(melted_spikeins.dt,
        endogenous_counts.df, counts2mpc.fit ,plot=plot,
        file=file, dropout_inflate=dropout_inflate,
        sample_number=sample_number)
    ##Prepare final output
    metadata(SCEList)<-list(dropout_parameters=noise_models$dropout_parameters,
        spikein_parameters=model_parameters, spikeConcentrations=spike_conc.df,
        genewiseDropouts=noise_models$genewise_dropouts)
    if (write.noise.model==TRUE) {
        write_noise_model(SCEList=SCEList, file=file)
    }
    SCEList
}
