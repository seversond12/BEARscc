# R script
# David T Severson
# Scope: Function in BEARscc used to estimate parameters and dropouts from ERCC spike-ins.

################################### INTERNAL FUNCTIONS ###############################################
##change to datatable, annotate with actual transcript count and optionally normalize by sequencing depth
prepare_data<-function(x,y){
    x$transcripts<-y[rownames(x),1]
    x
}

##computes linear fit for mean and standard deviation across ERCCs in the dataset
##x is prepared dataframe with cells as columns and ERCC ids as rows, each element is an observed unnormalized count
estimate_mu2sigma<-function(x, plot,sd_inflate, file="Rplot"){
  ##Correlate mean to standard deviation, so that negative binomial is a function of mean
  x1<-melt(data.table(x, keep.rownames = TRUE), id.vars = c("rn", "transcripts"))
  ERCC.sd<-data.table(x1[,sd(value),by=c("rn", "transcripts")])
  ERCC.sd$mean<-x1[,mean(value),by=c("rn", "transcripts")]$V1
  ERCC.sd[ERCC.sd$V1==0]$V1<-NA
  ERCC.sd[ERCC.sd$V1==0]$mean<-NA
  mu2sigma.fit<-lm(log2(ERCC.sd$V1)~log2(ERCC.sd$mean))
  ##Read correlation paramaters into data.frame
  mu2sigma<-data.frame(mu2sigma.slope=coef(mu2sigma.fit)[2], mu2sigma.intercept=coef(mu2sigma.fit)[1], max.res=max(mu2sigma.fit$residuals))
  if (plot==TRUE){
    g<-ggplot(data=ERCC.sd,aes(x=log2(mean), y=log2(V1)))+geom_point()+geom_abline(slope = coef(mu2sigma.fit)[2],intercept = coef(mu2sigma.fit)[1])+
      xlab("\nMean ERCC expression / log2")+ylab("Standard deviation of ERCC expression / log2\n")+
      ggtitle("Relationship between standard deviation and mean of ERCC spike in measurements\n")+
      geom_abline(slope=coef(mu2sigma.fit)[2],color="red", intercept=coef(mu2sigma.fit)[1]+sd_inflate*max(mu2sigma.fit$residuals))
    ggsave(paste(file,"mu2sigma.pdf", sep="_"), plot=g, device="pdf", width=10.5, height=7.5, units = "in")
  }
  mu2sigma
}

##function that computes linear fit parameters in order to model ERCC dispersion as proxy for technical noise in the dataset
##x is a dataframe with cells as columns and ERCC ids as rows, each element is an observed unnormalized count
##y is a dataframe with ERCC ids as rownames and the 1st column is the number of transcripts spiked into the sample
##granularity is the number of bins to divide the distributions into in order to determine maximal alpha
compute_alpha<-function(x,estimate_mu2sigma,sd_inflate,plot,granularity, file, alpha_granularity=alpha_granularity){
  mu2sigma<-estimate_mu2sigma(x, plot,sd_inflate, file=file)
  x1<-melt(data.table(x, keep.rownames = TRUE), id.vars = c("rn", "transcripts"))
 if (plot==TRUE){
   testing<-melt(data.table(x, keep.rownames = TRUE)[,-1, with=FALSE], id.vars = c("transcripts"))
    g<-ggplot(data=testing[transcripts>=1,], aes(x=factor(round(log10(transcripts+0.25), digits=3)),y=log2(value+0.25)))+geom_violin(scale = "width")+
      xlab("\nERCC molecules per cell / log10")+ylab("Measured counts / log2\n")+
      ggtitle("The count measurement distribution as a function of actual transcript number.\n")
    ggsave(paste(file,"countsdistribution_vs_moleculespercell.pdf", sep="_"), plot=g, device="pdf", width=10.5, height=7.5, units = "in")
  }
  y<-data.table(x1, key="rn")[,mean(value), by=c("rn", "transcripts")]
  setkey(x1, rn)
  colnames(y)<-c("rn","transcripts","Mean")
  y<-y[Mean>0,]
  ERCCnames<-as.character(y$rn)
  setkey(y, "rn")
  z2<-data.frame(as.character("DUMMY"), 1, 1, 1, 1)
  colnames(z2)<-c("rn","Mean", "p.nbinom","p.poisson","r.value")
  for (j in seq(from=0, to=1, by=alpha_granularity)){
    if (j %in% seq(from=0, to=1, by=0.05)){
      print(paste("Estimating error for ERCCs with alpha = ", j,sep=""))
    }
    k=1-j
    y3<-data.frame(as.character("DUMMY"), 1, 1, 1, 1)
    colnames(y3)<-c("rn","Mean", "p.nbinom","p.poisson","r.value")
    for (i in ERCCnames){
      x2<-x1[rn==i]
      p<-dim(x1[rn==i])[1]
      y.name<-as.character(y[rn==i]$rn)
      y.mean<-y[rn==i]$Mean
      hx1 <- rnbinom(p*j,size = mean(x2$value)^2/((2^(mu2sigma[1,1]*log2(mean(x2$value))+mu2sigma[1,2]+mu2sigma[1,3]))^2-mean(x2$value)),mu = mean(x2$value))
      hx2 <- rpois(p*k,mean(x2$value))
      if (j==0){
        hx.mixed<-hx2
      }else if (j==1){
        hx.mixed <- hx1
      }else {
        hx.mixed<-append(hx1, hx2)
      }
      brakes<-seq(from=0,to=max(hx.mixed, x2$value)+max(hx.mixed, x2$value)/granularity, by=max(hx.mixed, x2$value)/granularity)
      h1<-hist(x2$value, plot = FALSE, breaks =brakes)
      h.mixed<-hist(hx.mixed,plot = FALSE, breaks = brakes)
      rvalue<-sum(abs(h.mixed$counts-h1$counts)/(2*p))
      y2<-data.frame(rn=y.name,Mean=y.mean,p.nbinom=j,p.poisson=k,r.value=rvalue)
      y4<-rbind(y3,y2)
      y3<-y4
    }
    z1<-rbind(z2,y3)
    z2<-data.table(z1, key="rn")[!(rn=="DUMMY"),]
  }
  z2<-data.table(z2, key="rn")
  DT<-z2[,min(r.value), by=rn]
  y1<-data.frame(as.character("DUMMY"), 1, 1, 1, 1)
  colnames(y1)<-c("rn","Mean", "p.nbinom","p.poisson","r.value")
  for (i in DT$rn){
    y2<-z2[rn==i & r.value==DT[rn==i]$V1,]
    y3<-rbind(y1,y2)
    y1<-y3
  }
  rvalues.dt<-data.table(y1[-1], key="rn")
  rvalues.dt<-rvalues.dt[rvalues.dt$Mean>=1,]
  alpha2mean<-lm(rvalues.dt$p.nbinom~log2(rvalues.dt$Mean+0.025))
  if (plot==TRUE){
    g<-ggplot(data=rvalues.dt, aes(x=log2(Mean+0.025),y=p.nbinom))+geom_point(fill="red",color="black", size=2.2,alpha=0.6,  pch=21)+
      xlab("\nObserved mean ERCC expression, log2(expression)")+ylab("Neg. binomial contribution parameter, alpha\n")+
      ggtitle("The neg. binomial contribution is a function of gene expression.\n")+
      geom_abline(slope=coef(alpha2mean)[2], intercept = coef(alpha2mean)[1])
    ggsave(paste(file,"alpha2mu_correlation.pdf",sep="_"), plot=g, device="pdf", width=10.5, height=7.5, units = "in")
  }
  alpha2mu<-data.frame(alpha2mu.slope=coef(alpha2mean)[2], alpha2mu.intercept=coef(alpha2mean)[1], sd.inflate=sd_inflate)
  parameters<-data.frame(cbind(mu2sigma,alpha2mu), row.names=c("parameter.value"))
  parameters
}

##function that computes the nulls using the estimated parameters from compute_alpha funciton
##takes prepared dataframe of ERCCs as imput and parameters from compute_alpha
compute_models<-function(x, parameters){
  mu2sigma.slope<-parameters$mu2sigma.slope
  sd_inflate<-parameters$sd.inflate
  max.res<-parameters$max.res
  mu2sigma.intercept<-parameters$mu2sigma.intercept
  alpha2mu.slope<-parameters$alpha2mu.slope
  alpha2mu.intercept<-parameters$alpha2mu.intercept
  x2<-data.frame(rn="", transcripts=-1, variable="", value=-1)
  for (i in unique(x$rn)){
    x1<-x[rn==i]
    t<-x[rn==i]$transcripts
    y<-length(x1$value)
    p<-alpha2mu.slope*log2(mean(x1$value)+0.025)+alpha2mu.intercept
    if (p>1) {p=1}
    if (p<0) {p=0}
    hx.1 <- rnbinom(ceiling(y*p),size = mean(x1$value)^2/((2^(mu2sigma.slope*log2(mean(x1$value))+mu2sigma.intercept+max.res*sd_inflate))^2-mean(x1$value)),mu = mean(x1$value))
    hx.2 <- rpois(ceiling(y*(1-p)),mean(x1$value))
    hx1 <- rnbinom(y,size = mean(x1$value)^2/((2^(mu2sigma.slope*log2(mean(x1$value))+mu2sigma.intercept+max.res*sd_inflate))^2-mean(x1$value)),mu = mean(x1$value))
    hx2 <- rpois(y,mean(x1$value))
    hx<-sample(append(hx.1, hx.2),y,replace = FALSE)
    sims.m<-data.frame(rn=i,transcripts=t, variable="Mixed", value=hx)
    sims.nb<-data.frame(rn=i,transcripts=t, variable="NBinom", value=hx1)
    sims.p<-data.frame(rn=i,transcripts=t, variable="Poisson", value=hx2)
    x3<-rbind(x2,sims.m,sims.nb,sims.p, x1)
    x2<-x3
  }
  x2[-1,]
}

##function that computes the linear fit parameters in order to model the false zeroes
##x is dataframe with cells as columns and ERCC ids as rows, each element is an observed unnormalized count filtered for ERCCs with missed detection
estimate_undetected2molpercell<-function(x, plot, file){
  ##Correlate probability of missed detection to number of ERCC molecules per cell
  undetected<-data.table(x[,rowSums(.SD==0), .SD=c(2:(dim(x)[2]-1))]/(dim(x)[2]-2))
  undetected$rn<-as.character(x$rn)
  undetected$transcripts<-as.numeric(x$transcripts)
  undetected<-undetected[undetected$transcripts>0.9,]
  undetected2mpc.fit<-lm(undetected$V1~log2(undetected$transcripts))
  ##Read correlation paramaters into data.frame
  undetected2molpercell<-data.frame(undetected2mpc.slope=coef(undetected2mpc.fit)[2], undetected2mpc.intercept=coef(undetected2mpc.fit)[1])
  if (plot==TRUE){
    g<-ggplot(data=undetected,aes(x=log2(transcripts), y=V1))+geom_point()+
      xlab("\nERCC actual counts, log2(molecules per cell)")+ylab("Fraction of cells where observed ERCC count is zero")+
      ggtitle("Conditional probability of observing a zero given transcript concentration\n")+
      geom_abline(slope=coef(undetected2mpc.fit)[2],color="red", intercept=coef(undetected2mpc.fit)[1])
    ggsave(paste(file,"undetected2molpercell.pdf", sep="_"), plot=g, device="pdf", width=11, height=8.5, units = "in")
  }
  undetected2molpercell
}

##function that computes the linear fit parameters from count data to actual molecules per cell
counts2mpc<-function(x,plot, file){
  x<-melt(data.table(x, keep.rownames = TRUE)[,-1, with=FALSE], id.vars = c("transcripts"))
  x<-x[transcripts>0.9,][,mean(value),by=transcripts]
  counts2mpc.fit<-lm(log2(x$V1+0.025)~log2(x$transcripts))
  if (plot==TRUE){
    g<-ggplot(data=x, aes(x=log2(transcripts),y=log2(V1+0.025)))+geom_point()+
      xlab("\nERCC molecules per cell, log2(molecules)")+ylab("Observed counts, log2(counts) \n")+
      ggtitle("Correlation of measured counts as a function of actual transcript number.\n")+
      geom_abline(slope=coef(counts2mpc.fit)[2], intercept=coef(counts2mpc.fit)[1], color="red")
    ggsave(paste(file,"cor_counts2moleculespercell.pdf", sep="_"), plot=g, device="pdf", width=10.5, height=7.5, units = "in")
  }
  counts2mpc.fit
}

##function that computes genewise probabilities that observed values are zero
compute_genewise_zeroinflation<-function(x){
  bayes_list<-list()
  bayes_list$genewise<-data.frame(rowSums(x==0)/dim(x)[2])
  bayes_list
}

### By Bayes law, P(x=k | y= 0) = (P(y=0 | x=k) * P(y=0))/P(x=k), where k=1...34 in N
### Each transcript is only present in the set once, so the distribution
###for K random variable is uniform, where P(x=k)= 1/34 = 0.02941176
### By the law of total probability, P(y=0) = sum,from k=1 to 34,{P(y=0|x=k)*P(x=ki)}
## This is equivalent to sum,from k=1 to 34,{P(y=0|x=k)} * 1/34,
### x is the prepared dataframe of ERCC spikeins, y is the dataframe of raw counts of actual cellular expression
estimate_missingdata<-function(x, y, counts2mpc.fit, plot, file, dropout_inflate=dropout_inflate){
  bayes_list<-compute_genewise_zeroinflation(y)
  ERCC.dt<-data.table(x, keep.rownames = TRUE)
  ERCC.dt<-ERCC.dt[ERCC.dt[,rowSums(.SD==0)>1, .SD=c(2:dim(ERCC.dt)[2])],][transcripts>0,]
  undetected2mpc<-estimate_undetected2molpercell(ERCC.dt, plot=plot, file=file)
  kmax<-ceiling(2^(-as.numeric(undetected2mpc[2])/((1/as.numeric(dropout_inflate))*as.numeric(undetected2mpc[1]))))
  k<-seq(from=0, to=kmax, by=1)
  y0.df<-data.frame(k=k)
  y0.df$Py0givenk<-(1/as.numeric(dropout_inflate))*as.numeric(undetected2mpc[1])*log2(y0.df$k)+as.numeric(undetected2mpc[2])
  y0.dt<-data.table(y0.df, keep.rownames = TRUE)
  y0.dt[y0.dt$Py0givenk>1 | y0.dt$Py0givenk==Inf,]$Py0givenk<- 1 ## Assume that if k=0, then the probability of observing 0 is 1
  y0.dt[y0.dt$Py0givenk<0,]$Py0givenk<-0 ## squeeze negative to 0
  y0.dt<-y0.dt[y0.dt$Py0givenk>0,]
  y0.dt$Counts<-2^(coef(counts2mpc.fit)[2]*log2(y0.dt$k)+coef(counts2mpc.fit)[1])
  size<-dim(y0.dt)[1]-1
  colnames(bayes_list$genewise)<-"py0"
  bayes_list$genewise$pk0<-(bayes_list$genewise$py0-sum(y0.dt$Py0givenk[-1])/size)/(1-sum(y0.dt$Py0givenk[-1])/size)
  if (length(bayes_list$genewise[bayes_list$genewise$pk0<0,]$pk0)>0) {
    bayes_list$genewise[bayes_list$genewise$pk0<0,]$pk0<-0
    bayes_list$genewise[bayes_list$genewise$pk0==0 & !(bayes_list$genewise$py0==0),]$py0<-sum(y0.dt$Py0givenk[-1])/size
  }
  bayes_list$genewise$pki<-(1-bayes_list$genewise$pk0)*1/size
  py0givenk<-as.vector(y0.dt$Py0givenk[-1])
  gene.pki<-as.vector(bayes_list$genewise$pki)
  names(gene.pki)<-as.character(rownames(bayes_list$genewise))
  bayes_list$inferred_prob<-py0givenk%*%t(gene.pki)
  bayes_list$inferred_prob<-sweep(rbind(t(bayes_list$genewise$pk0),bayes_list$inferred_prob), 2, bayes_list$genewise$py0, "/")
  bayes_list$inferred_prob<-data.table(data.frame(bayes_list$inferred_prob), keep.rownames = TRUE)
  bayes_list$bayes_parameters<-y0.dt
  setkey(bayes_list$bayes_parameters,rn)
  setkey(bayes_list$inferred_prob, rn)
  bayes_list$bayes_parameters<-bayes_list$bayes_parameters[bayes_list$inferred_prob,]
  setkey(bayes_list$bayes_parameters, k)
  bayes_list
}

################################## USER COMMANDS ############################################
##x is the dataframe of ERCCs, y is the data frame of ERCCs with actual molecules per cell in the second column, z is the counts matrix
##model_view should be a vector of c("Optimized", "Poisson", "Neg. Binomial", "Observed")
estimate_noiseparameters<-function(spike_counts.df,endogenous_counts.df,spike_conc.df, plot=FALSE, sd_inflate=0,granularity=300,write.noise.model=TRUE,file="noise_estimation",model_view=c("Observed","Optimized"),total_sampling=2500, dropout_inflate=1, alpha_granularity=0.005){
  force(granularity)
  force(sd_inflate)
  force(alpha_granularity)
  ERCC.prepared.counts<-prepare_data(spike_counts.df, spike_conc.df)
  actual.prepared.counts<-endogenous_counts.df
  print("Fitting parameter alpha to establish ERCC-derived noise model.")
  parameters<-compute_alpha(ERCC.prepared.counts,estimate_mu2sigma, granularity=granularity, sd_inflate=sd_inflate, plot=plot, file=file, alpha_granularity=alpha_granularity)
  ERCC.m.counts<-melt(data.table(ERCC.prepared.counts, keep.rownames = TRUE), id.vars = c("rn", "transcripts"))
  ERCC.m.valid<-ERCC.m.counts[,mean(value)>0, by="rn"]
  setkey(ERCC.m.valid, "rn")
  setkey(ERCC.m.counts,"rn")
  ERCC.m.counts<-ERCC.m.valid[ERCC.m.counts,][V1==TRUE,][,-2,with=FALSE]
  models.dt<-compute_models(ERCC.m.counts, parameters)
  models.dt$distribution<-gsub(TRUE,"Observed",models.dt$variable != "Mixed" & models.dt$variable != "NBinom" & models.dt$variable !="Poisson")
  models.dt[models.dt$variable=="Mixed",]$distribution<-"Optimized"
  models.dt[models.dt$variable=="Poisson"]$distribution<-"Poisson"
  models.dt[models.dt$variable=="NBinom"]$distribution<-"Neg. Binomial"
  if (plot==TRUE){
    ERCCnames<-as.character(unique(models.dt[order(transcripts)]$rn))
    round(length(ERCCnames)/6)
    dir.create("ERCCfit_plots")
    for (i in seq(from=1, to=round(length(ERCCnames)/6), by=1)){
      p<-ERCCnames[seq(from=i, to=length(ERCCnames), by=round(length(ERCCnames)/6))]
      g<-ggplot(data=models.dt[rn %in% p,][distribution %in% model_view], aes(x=value))+geom_histogram(aes(fill=factor(distribution, levels=c("Poisson", "Neg. Binomial", "Optimized", "Observed"))),alpha=0.30,position="identity", bins=granularity)+
        facet_wrap(~rn, ncol=3,nrow=4, scales="free")+ylab("Observation density\n")+labs(fill="Model")+xlab("\nERCC expression / normalized counts")+ggtitle("Approximating technical noise using a refined mixed distribution.\n")+scale_color_discrete(guide=FALSE)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10))
      savefile<-paste("ERCCfit_plots/",paste(file,"ERCCfits_subset", i, "plot.pdf",sep="_"),sep="")
      ggsave(savefile, plot=g, device="pdf", width=10.5, height=7.5, units = "in")
    }
  }
  counts2mpc.fit<-counts2mpc(ERCC.prepared.counts, plot=plot, file=file)
  bayes_dropouts<-estimate_missingdata(ERCC.prepared.counts, actual.prepared.counts, counts2mpc.fit ,plot=plot, file=file, dropout_inflate=dropout_inflate)
  bayes_dropouts$ERCC_parameters<-parameters
  noise_model<-bayes_dropouts
  noise_model$models.dt<-models.dt
  noise_model$original.counts<-endogenous_counts.df
  noise_model$spike.conc<-spike_conc.df
  noise_model$spike.counts<-spike_counts.df
  if (write.noise.model==TRUE) {
    write.table(data.table(noise_model$ERCC_parameters, keep.rownames = TRUE),file=paste(file, "parameters4randomize.xls", sep="_"), quote=FALSE, sep="\t", row.names = FALSE)
    write.table(noise_model$bayes_parameters,file=paste(file, "bayesianestimates.xls", sep="_"), quote=FALSE, sep="\t", row.names=FALSE )
    print("Parameters have been saved into the current directory in bayesianestimates.xls and parameters4randomize.xls.")
    noise_model
  } else {
    noise_model
  }
}
