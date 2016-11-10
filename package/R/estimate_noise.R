# R script
# David Severson
# Scope: To estimate parameters and generate noisy counts tables. Package allows users to vet clusters

#### internal functions##################################################
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
compute_alpha<-function(x,estimate_mu2sigma,sd_inflate,plot,granularity, file){
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
  for (j in seq(from=0, to=1, by=0.005)){
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
  print("Alpha with least error in fit has been computed for each measured ERCC.")
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
  print("Model parameters have been estimated. They are: ")
  print(parameters)
  parameters
}

##function that computes the nulls using the estimated parameters from compute_alpha funciton
##takes prepared dataframe of ERCCs as imput and parameters from compute_alpha
compute_models<-function(x, parameters){
  print("We are now sampling from distributions to create vector for plotting")
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
  print("Determing ERCCs with dropout observations in actual detection")
  print(paste("with dropout inflation, ", dropout_inflate, ".", sep=""))
  bayes_list<-compute_genewise_zeroinflation(y)
  ERCC.dt<-data.table(x, keep.rownames = TRUE)
  ERCC.dt<-ERCC.dt[ERCC.dt[,rowSums(.SD==0)>1, .SD=c(2:dim(ERCC.dt)[2])],][transcripts>0,]
  undetected2mpc<-estimate_undetected2molpercell(ERCC.dt, plot=plot, file=file)
  kmax<-ceiling(2^(-as.numeric(undetected2mpc[2])/((1/as.numeric(dropout_inflate))*as.numeric(undetected2mpc[1]))))
  print(kmax)
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


##actually applies the estimated noise to the counts table
randomizer<-function(x,parameters,total_sampling){
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
permute_count<-function(x, probs4detection.k, probabilityA, parameters, gene, randomizer, total_sampling){
  force(parameters)
  if (x==0) {
    nx<-sample(as.numeric(rownames(probs4detection.k)),1,prob=probabilityA, replace=TRUE)
    probabilityB<-probs4detection.k[as.character(ceiling(nx)),]$Py0givenk
    nx<-sample(c(0,nx),1, prob=as.numeric(c(probabilityB, 1-probabilityB)))
    if (nx>0) {
      nx<-randomizer(probs4detection.k[as.character(nx),]$Counts, parameters, total_sampling)
    }
  }
  else {
    if (ceiling(x) < max(as.numeric(rownames(probs4detection.k)))){
      probabilityB<-as.numeric(probs4detection.k[as.character(ceiling(x)),]$Py0givenk)
      nx<-sample(c(0,x),1, prob=as.numeric(c(probabilityB, 1-probabilityB)))
      if (nx>0) {
        nx<-randomizer(nx, parameters, total_sampling)
      }
    } else {
      nx<-randomizer(x, parameters, total_sampling)
    }
  }
  nx
}

##seperates out genes for analysis
genewise_permute_count<-function(x,probs4detection.k, probs4detection.genes,parameters,randomizer,total_sampling=2500){
  probabilityA<-probs4detection.genes[gsub("-",".",x[1]),]
  force(total_sampling)
  apply(data.frame(as.numeric(x[-1])),1, `permute_count`, probs4detection.k,probabilityA, parameters,gene=gene_name,randomizer,total_sampling)
}

##create cell association matrices for compute consensus
create_cm<-function(z, names){
  names(z)<-names
  sapply(z, function(x){sapply(z, function(y){as.numeric(x==y)})})
}

###underlying function calculating cluster stability

##functions for cluster stability from estimate_noise.R
###x is the consensus matrix
###y is a dataframe of cluster identities for various k

calculate_cluster_metrics_by_cluster<-function(i,y,x){
  in.cluster<-as.character(names(y[y==i]))
  out.cluster<-as.character(names(y[!(y==i)]))
  if (length(in.cluster)>1){
    inside.stability<-(sum(x[in.cluster,in.cluster, drop=FALSE])-length(in.cluster))/(length(in.cluster)^2-length(in.cluster))
  } else {
    inside.stability<-1
  }
  promiscuity.matrix<-as.matrix(x[in.cluster,out.cluster, drop=FALSE])
  if (length(in.cluster) >= length(out.cluster)){
    outside.promiscuity<-mean(promiscuity.matrix)
  } else {
    prom.mean.dt<-data.table(data.frame(colMeans(promiscuity.matrix)), keep.rownames = TRUE)
    colnames(prom.mean.dt)[2]<-"mean.prom"
    prom.mean.dt<-prom.mean.dt[order(-mean.prom)]
    prom.mean.dt<-prom.mean.dt[1:length(in.cluster),]
    outside.promiscuity<-mean(promiscuity.matrix[,as.character(prom.mean.dt$rn), drop=FALSE])
  }
  score<-inside.stability-outside.promiscuity
  metrics.dt<-data.table(data.frame(score=score, Promiscuity=outside.promiscuity, Stability=inside.stability, size=length(in.cluster)))
  metrics.dt
}

calculate_cluster_metrics<-function(y,x,z){
  names(y)<-z
  cluster_names<-as.vector(as.character(unique(y)))
  names(cluster_names)<-as.character(unique(y))
  metric_list<-lapply(cluster_names, `calculate_cluster_metrics_by_cluster`, x=x, y=y)
  metric_list<-melt(metric_list, id.vars=c("size"))
  colnames(metric_list)[c(2,4)]<-c("metric","rn")
  metric_list
}

calculate_cell_metrics_by_cluster<-function(i, y, x){
  in.cluster<-as.character(names(y[y==i]))
  out.cluster<-as.character(names(y[!(y==i)]))
  if (length(in.cluster) < length(out.cluster)){
    if (length(in.cluster)==1){
      inside.stability<-1
      outside.promiscuity<-max(x[in.cluster, out.cluster])
      score<-1-outside.promiscuity
      cluster_metric.dt<-data.frame(cell=as.character(in.cluster[1]),score=score, Promiscuity=outside.promiscuity, Stability=inside.stability, cluster.size=length(in.cluster))
      cluster_metric.dt$cluster.size<-length(in.cluster)
    } else {
      stability.dt<-data.table(data.frame(x[in.cluster,in.cluster, drop=FALSE]), keep.rownames = TRUE)
      stability.dt<-stability.dt[,.(Stability = (rowSums(.SD)-1)/(length(in.cluster)-1)), keyby=rn]
      promiscuity.dt<-data.table(data.frame(x[in.cluster,out.cluster, drop=FALSE]), keep.rownames = TRUE)
      promiscuity.dt<-promiscuity.dt[,.(Promiscuity = apply(.SD,1, function(x, y=length(in.cluster)) {mean(sort(x, decreasing = TRUE)[1:y])})), keyby=rn]
      cluster_metric.dt<-stability.dt[promiscuity.dt]
      colnames(cluster_metric.dt)[1]<-"cell"
      cluster_metric.dt$score<-cluster_metric.dt$Stability-cluster_metric.dt$Promiscuity
      cluster_metric.dt$cluster.size<-length(in.cluster)
    }
  } else if (length(in.cluster) >= length(out.cluster)) {
    if (length(in.cluster)==dim(x)[1]) {
      stability.dt<-data.table(data.frame(x[in.cluster,in.cluster, drop=FALSE]), keep.rownames = TRUE)
      stability.dt<-stability.dt[,.(Stability = (rowSums(.SD)-1)/(length(in.cluster)-1)), by=rn]
      stability.dt$Promiscuity<-0
      colnames(stability.dt)[1]<-"cell"
      cluster_metric.dt<-stability.dt
      cluster_metric.dt$cluster.size<-length(in.cluster)
      cluster_metric.dt$score<-cluster_metric.dt$Stability-cluster_metric.dt$Promiscuity
    } else {
      stability.dt<-data.table(data.frame(x[in.cluster,in.cluster,drop=FALSE]), keep.rownames = TRUE)
      stability.dt<-stability.dt[,.(Stability = (rowSums(.SD)-1)/(length(in.cluster)-1)), keyby=rn]
      promiscuity.dt<-data.table(data.frame(x[in.cluster,out.cluster, drop=FALSE]), keep.rownames = TRUE)
      promiscuity.dt<-promiscuity.dt[,.(Promiscuity = rowMeans(.SD)), keyby=rn]
      cluster_metric.dt<-stability.dt[promiscuity.dt]
      cluster_metric.dt$score<-cluster_metric.dt$Stability-cluster_metric.dt$Promiscuity
      cluster_metric.dt$cluster.size<-length(in.cluster)
      colnames(cluster_metric.dt)<-c("cell", "Stability", "Promiscuity", "score", "cluster.size")
    }
  }
  cluster_metric.dt
}

calculate_cell_metrics<-function(y,x,z){
  names(y)<-z
  cluster_names<-as.vector(as.character(unique(y)))
  names(cluster_names)<-as.character(unique(y))
  cell_metric_list<-data.table(melt(lapply(cluster_names, `calculate_cell_metrics_by_cluster`, x=x, y=y), id.vars=c("cell", "cluster.size")))
  colnames(cell_metric_list)[c(3,5)]<-c("metric","rn")
  cell_metric_list
}

################################################################################

## user commands ###############################################################
##x is the dataframe of ERCCs, y is the data frame of ERCCs with actual molecules per cell in the second column, z is the counts matrix
##model_view should be a vector of c("Optimized", "Poisson", "Neg. Binomial", "Observed")
estimate_noiseparameters<-function(x,y,z, plot=FALSE, sd_inflate=0, normseq=TRUE,granularity=300,write.noise.model=TRUE,file="noise_estimation",model_view=c("Observed","Optimized"),total_sampling=2500, dropout_inflate=1){
  force(granularity)
  force(sd_inflate)
  print("Preparing data for parameter estimation.")
  ERCC.prepared.counts<-prepare_data(x, y)
  actual.prepared.counts<-z
  print("Fitting the linear model as a function of gene expression from ERCCs.")
  parameters<-compute_alpha(ERCC.prepared.counts,estimate_mu2sigma, granularity=granularity, sd_inflate=sd_inflate, plot=plot, file=file)
  if (plot==TRUE){
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
  print("Calculating the fit between actual transcripts and the counts observed.")
  counts2mpc.fit<-counts2mpc(ERCC.prepared.counts, plot=plot, file=file)
  bayes_dropouts<-estimate_missingdata(ERCC.prepared.counts, actual.prepared.counts, counts2mpc.fit ,plot=plot, file=file, dropout_inflate=dropout_inflate)
  bayes_dropouts$ERCC_parameters<-parameters
  noise_model<-bayes_dropouts
  noise_model$models.dt<-models.dt
  if (write.noise.model==TRUE) {
    write.table(data.table(noise_model$ERCC_parameters, keep.rownames = TRUE),file=paste(file, "parameters4randomize.xls", sep="_"), quote=FALSE, sep="\t", row.names = FALSE)
    write.table(noise_model$bayes_parameters,file=paste(file, "bayesianestimates.xls", sep="_"), quote=FALSE, sep="\t", row.names=FALSE )
    print("Parameters have been saved into the current directory in bayesianestimates.xls and parameters4randomize.xls.")
    noise_model
  } else {
    "Estimations are being printed"
    noise_model
  }
}

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

##executes noise perturbation with estimated parameters
create_noiseinjected_counts<-function(x,probs4detection, parameters,genewise_permute_count=genewise_permute_count,permute_count=permute_count, randomizer=randomizer,total_sampling=2500){
  probs4detection.genes<-t(data.frame(probs4detection, row.names = "k")[,4:eval(dim(probs4detection)[2]-1)])
  probs4detection.k<-data.frame(probs4detection[,2:4, with=FALSE],row.names = "k")
  noisy_counts<-data.table(x, keep.rownames = TRUE)[,apply(.SD,1 ,`genewise_permute_count`,probs4detection.k=probs4detection.k,probs4detection.genes=probs4detection.genes,total_sampling, randomizer=randomizer, parameters=parameters)]
  noisy_counts
}

##compute consensus from dataframe of perturbaton reclustering labels
compute_consensus<-function(x, plot=FALSE, file="Rplot", orig.labels=NA, num.clust=5){
  print("Calculating consensus matrix...")
  cm<-lapply(x, `create_cm`,names=rownames(x))
  cm<-Reduce("+",cm)/dim(x)[2]
  print("Calculating noise consensus clusterings...")
  rclust<-hclust(as.dist(abs(1-cm)), method="average")
  cclust<-hclust(as.dist(abs(1-t(cm))), method="average")
  clusts <- cutree(cclust, k = num.clust)
  clusts<-clusts[colnames(cm)]
  nccolors<-c("darkred", "dodgerblue", "blue3", "limegreen", "darkorange4", "darkcyan", "purple", "gold3","#dad1bf", "#5e4400", "#00abff", "#73a3df", "#d6d1c8","#514627", "maroon3", "#87d6ff", "#004eae", "#a4a097", "#efcf72", "#93d6ff", "#304865", "#76726a", "#b99e4e", "#0077c6", "#7da3d6", "#4a463f", "#847035", "#0077cf", "#55749c", "#e0d0ac", "#52461f", "#0f76be", "#84a2cd", "#ddd0b5", "#f2ce5c", "#004a7d", "#0050bf", "#ab9f86", "#ba9e43", "#9dd5ff", "#5d7393", "#a8a08f", "#85702c", "#31a6ff", "#38485e", "#797161", "#544516", "#004a85", "#b6d4ff", "#e6d099", "#f5ce44", "#0078d8", "#0052d0", "#e3d0a3", "#bd9e2b", "#0079e1", "#8ba2c4", "#ae9f7d", "#867020", "#004b8d", "#63738b", "#7c7159", "#8a7000", "#154975", "#3e4756", "#4d4637", "#fbce00", "#a6d5ff", "#bcd3f7", "#ebcf86", "#bf9d00", "#007ae9", "#c2d3ed", "#e9cf8f", "#564500", "#2f76b5", "#92a2bb", "#b39f6a", "#d1d1d1", "#004c95", "#697383", "#b19f74", "#a0a0a0", "#007cff", "#c8d2e4", "#807148", "#727272", "#5ba5f1", "#0054e1", "#7e7150", "#474747", "#4075ad", "#97a1b2", "#4f462f", "#00ddff", "#004c9d", "#9ca1a9", "#edcf7c", "#29d9ff", "#26496d", "#6e727a", "#b79e58", "#50d8ff", "#aed4ff", "#43474e", "#b59e61", "#67d7ff", "#69a4e8", "#cdd2db", "#82703f", "#79d7ff", "#4c74a4")
  newclusters<-as.character(clusts[colnames(cm)])
  if (plot==TRUE && is.na(orig.labels)){
    print("Creating plot without original cluster labels, because none were provided.")
    aheatmap(cm,treeheight=80, annCol = list("Consensus Clusters"=newclusters),labCol=NA, labRow=" ", Rowv = rclust, Colv=cclust)
  } else if (plot==TRUE && !(is.na(orig.labels))) {
    print("Creating plot with original cluster labels, because they were provided")
    oldclusters<-as.character(orig.labels[colnames(cm),])
    aheatmap(cm,treeheight=80, annColors =list("Consensus Clusters"=rev(nccolors), "Original Clusters"=nccolors) ,annCol = list("Consensus Clusters"=newclusters, "Original Clusters"=oldclusters), labCol=NA, labRow=" ", Rowv = rclust, Colv=cclust)
  }
  clusts<-data.table(data.frame(clusts), keep.rownames = TRUE)[order(clusts)]
  list(consensus_matrix=cm,clusts=clusts)
}


cluster_consensus<-function(x,y, method="complete"){
  cclust<-hclust(as.dist(abs(1-t(x))), method=method)
  clusts<-data.frame(cutree(cclust, k=y))
  clusts
}


##seperately plot cell and cluster noise robustness index
# y is a data frame of the final clustering results perturbed labels
# x is a data frame of the consensus matrix. need to make robust to bug when names start with numbers...
report_cell_metrics<-function(y,x){
  z<-rownames(x)
  y$X1<-rep("1",dim(y)[1])
  y<-data.table(y)
  cell.scores<-list()
  cell.scores<-lapply(y, `calculate_cell_metrics`,x=x, z=z)
  cell.scores<-data.table(melt(cell.scores, id.var=c("rn","cell","cluster.size", "metric")))
  cell.scores$L1<-as.factor(gsub("^X","", cell.scores$L1))
  cell.scores<-cell.scores[,-5, with=FALSE]
  cell.scores
}

###x is the consensus matrix with rows and columns as cell names
###y is a dataframe of cluster identities for various k in columns and rownames as cells  
report_cluster_metrics<-function(y,x, weighted_mean=FALSE,plot=FALSE, file="Rplot"){
  print(weighted_mean)
  z<-rownames(y)
  y<-data.table(y)
  cluster.scores<-list()
  oneclust.stability<-(sum(x)-dim(x)[1])/dim(x)[1]^2
  cluster.scores$X1<-melt(data.table(rn="1",score=oneclust.stability,Promiscuity=0,Stability=oneclust.stability, size=dim(x)[1]), id.vars=c("rn", "size"))
  colnames(cluster.scores$X1)[3]<-"metric"
  cluster.scores<-append(cluster.scores,lapply(y, `calculate_cluster_metrics`,x=x, z=z))
  cluster.scores.m<-data.table(melt(cluster.scores, id.var=c("rn","size","metric")))
  cluster.scores.m<-cluster.scores.m[,-4, with=FALSE]
  cluster.scores.m$L1<-gsub("^X","", cluster.scores.m$L1)
  cluster.scores.m$singlet<-cluster.scores.m$size==1
  cluster.scores.m$singlet<-factor(gsub(FALSE,"Cell number > 1", gsub(TRUE, "Cell number = 1", cluster.scores.m$singlet)),levels = c("Cell number > 1","Cell number = 1"))
  if (weighted_mean==TRUE){
    mean.overall<-cluster.scores.m[,sum(value*size/as.numeric(cluster.scores.m[L1=="1",]$size[1])), by=c("L1", "metric")]
  }else {
    mean.overall<-cluster.scores.m[,mean(value), by=c("L1", "metric")]
  }
  colnames(mean.overall)[3]<-"Overall.mean"
  setkeyv(cluster.scores.m,c("L1", "metric"))
  setkeyv(mean.overall, c("L1", "metric"))
  cluster.scores.m<-cluster.scores.m[mean.overall]
  numbers<-as.character(sort(unique(na.omit(as.numeric(gsub("^X","",cluster.scores.m$L1))))))
  other<-unique(cluster.scores.m$L1[!(gsub("X","",cluster.scores.m$L1)%in%numbers)])
  level<-as.character(c(numbers,other))
  cluster.scores.m$L1<-factor(cluster.scores.m$L1, levels = level)
  if (plot==TRUE){
    g<-ggplot(data=cluster.scores.m, aes(x=L1, y=value))+geom_boxplot(color="grey20", fill="grey85")+
      guides(fill="none", color="none")+xlab("\nCluster number")+ylab("Cluster metric\n")+facet_wrap(~metric, ncol=1,nrow=3, scales = "free")+
      ggtitle("Cluster technical consensus scores for various clustering\n")+labs(alpha="Cluster size")+
      geom_point(aes(x=factor(L1), y=Overall.mean),fill="red", size=2.5, color="black", pch=23)+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
    savefile<-paste(file, "_cluster_scores.pdf", sep="")
    ggsave(savefile, plot=g, device="pdf", width=8.5, height=11, units = "in") 
  }
  cluster.scores.m
}
