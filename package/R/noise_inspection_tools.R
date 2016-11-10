# R script
# David T Severson
# Scope: Functions in BEARscc used to visualize noise estimated from ERCCs in single cell data.

###################### USER FUNCTIONS ##############################
##seperately plot deviation of fitted model (or standard poisson/neg. binomial from observed counts)
##needs editing....
visualize_model_error<-function(x,y,granularity=300,parameters){
  fit.slope.mean2sd<-as.numeric(parameters[1])
  fit.intercept.mean2sd<-as.numeric(parameters[2])
  fit.slope.p2mean<-as.numeric(parameters[2])
  sd_inflate<-as.numeric(parameters[6])
  max.res=max(SD2mean.fit$residuals)
  brakes<-seq(from=0,to=round(max(x$value)*1.2), by=1300/granularity)
  y3<-data.frame(as.character("DUMMY"), 1, 1, 1,1)
  colnames(y3)<-c("rn","Mean", "nbinom.sd", "poisson.sd", "mixed.sd")
  for (i in ERCCnames){
    x1<-x[rn==i]
    y.name<-as.character(y[rn==i]$rn)
    y.mean<-y[rn==i]$Mean
    hx1 <- rnbinom(1107,size = mean(x1$value)^2/((2^(fit.slope.mean2sd*log2(mean(x1$value))+fit.intercept.mean2sd))^2-mean(x1$value)),mu = mean(x1$value))
    hx2 <- rpois(1107,mean(x1$value))
    p<-1-fit.slope.p2mean*y.mean
    if (p>1) {p=1}
    if (p<0) {p=0}
    hx.1 <- rnbinom(1107*p,size = mean(x1$value)^2/((2^(fit.slope.mean2sd*log2(mean(x1$value))+fit.intercept.mean2sd+sd_inflate*max.res))^2-mean(x1$value)),mu = mean(x1$value))
    hx.2 <- rpois(1107*(1-p),mean(x1$value))
    h.mixed<-append(hx.1,hx.2)
    h.mixed<-hist(h.mixed, plot=FALSE, breaks=brakes)
    h1<-hist(x1$value,plot = FALSE, breaks=brakes)
    h.nbinom<-hist(hx1,plot = FALSE, breaks=brakes)
    h.poisson<-hist(hx2,plot = FALSE, breaks=brakes)
    sd.nbinom<-sum(abs(h.nbinom$counts-h1$counts)/2214)
    sd.poisson<-sum(abs(h.poisson$counts-h1$counts)/2214)
    sd.mixed<-sum(abs(h.mixed$counts-h1$counts)/2214)
    y2<-data.frame(rn=y.name,Mean=y.mean,nbinom.sd=sd.nbinom,poisson.sd=sd.poisson, mixed.sd=sd.mixed)
    y4<-rbind(y3,y2)
    y3<-y4
  }
  rvalues.dt<-melt(data.table(y3[-1,], key="rn"), id.vars=c("rn", "Mean"))
  rvalues.dt<-compute_fractionvalues(ERCC.counts,ERCC.means.dt,50)
  ggplot(data=rvalues.dt, aes(x=log2(Mean),y=value))+geom_point(aes(fill=gsub("nbinom.sd","Neg. Binomial",gsub("poisson.sd","Poisson",gsub("mixed.sd","Mixed",variable)))),color="black", size=2.2,alpha=0.6,  pch=21)+
    xlab("\nObserved mean ERCC expression / log2")+ylab("Average deviation\n")+labs(fill="Random variable")+
    ggtitle("The goodness of the statistical fit for each model as a function gene expression.\n")

}

##seperately plot observed spread of ERCCs, could be useful to compare normalized and unnormalized dispersion
##needs editing...
plot_obs_spread<-function(x,y, file="Rplot",parameters){
  ERCC.prepared.counts<-prepare_data(x, y)
  ERCC.m.counts<-melt(data.table(ERCC.prepared.counts, keep.rownames = TRUE), id.vars = c("rn", "transcripts"))
  ERCC.m.valid<-ERCC.m.counts[,mean(value)>0, by="rn"]
  setkey(ERCC.m.valid, "rn")
  setkey(ERCC.m.counts,"rn")
  ERCC.m.counts<-ERCC.m.valid[ERCC.m.counts,][V1==TRUE,][,-2,with=FALSE]
  models.dt<-compute_models(ERCC.m.counts, parameters)
  models.dt$distribution<-gsub(TRUE,"Observed",models.dt$variable != "Mixed" & models.dt$variable != "NBinom" & models.dt$variable !="Mixed")
  models.dt[models.dt$variable=="Mixed",]$distribution<-"Optimized"
  models.dt[models.dt$variable=="Poisson"]$distribution<-"Poisson"
  models.dt[models.dt$variable=="NBinom"]$distribution<-"Neg. Binomial"
  ERCCnames<-as.character(unique(models.dt[order(transcripts)]$rn))
  round(length(ERCCnames)/6)
  dir.create("ERCCfit_plots")
  for (i in seq(from=1, to=round(length(ERCCnames)/6), by=1)){
    p<-ERCCnames[seq(from=i, to=length(ERCCnames), by=round(length(ERCCnames)/6))]
    g<-ggplot(data=models.dt[rn %in% p,][distribution %in% model_view], aes(x=value))+geom_histogram(aes(color=factor(distribution, levels=c("Poisson", "Neg. Binomial", "Optimized", "Observed")), fill=factor(distribution, levels=c("Poisson", "Neg. Binomial", "Optimized", "Observed"))), stat="identity", alpha=0.2)+
      facet_wrap(~rn, ncol=3,nrow=4, scales="free")+ylab("Density\n")+labs(fill="Data source")+xlab("\nERCC expression / normalized counts")+ggtitle("Approximating technical noise using a refined mixed distribution.\n")+scale_color_discrete(guide=FALSE)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10))
    savefile<-paste("ERCCfit_plots/",paste(file,"ERCCfits_subset", i, "plot.pdf",sep="_"),sep="")
    ggsave(savefile, plot=g, device="pdf", width=10.5, height=7.5, units = "in")
  }
}

##requires x bayesian list, and a gene list
##needs editing
plot_bayes_inference<-function(x){
  y0.dt$sep<-gsub("TRUE","k=0",gsub("FALSE","k>0",y0.dt$pkgiveny0>0.1))
  g<-ggplot(data=y0.dt, aes(x=k,y=pkgiveny0))+geom_point() +labs(fill="Detected?")+
    xlab("\nk, theoretical number of ERCC transcripts in sample")+ylab("P(theor. counts = k | obs. counts = 0)\n")+
    ggtitle("Probability that the theoretical number of transcripts is k, given the observed count value is 0.")+
    facet_wrap(~sep, scales="free_y", nrow = 2,ncol = 1, as.table=FALSE)
  ggsave(paste(file,"dropoutrate_vs_transcriptsnumber.pdf", sep="_"), plot=g, device="pdf", width=10.5, height=7.5, units = "in")
  g<-ggplot(data=y0.dt, aes(x=Counts,y=pkgiveny0))+geom_point() +labs(fill="Detected?")+
    xlab("\nCounts, number of ERCC transcripts in sample")+ylab("P(true counts = k | obs. counts = 0)\n")+
    ggtitle("Probability that the observed count value should have been k, given the observed count value is 0.")+
    facet_wrap(~sep, scales="free_y", nrow = 2,ncol = 1, as.table=FALSE)
  ggsave(paste(file,"dropoutrate_vs_count.pdf", sep="_"), plot=g, device="pdf", width=10.5, height=7.5, units = "in")
}
