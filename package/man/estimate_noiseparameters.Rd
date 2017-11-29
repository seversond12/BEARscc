\name{estimate_noiseparameters}
\alias{estimate_noiseparameters}
\alias{prepare_data}
\alias{estimate_mu2sigma}
\alias{compute_alpha}
\alias{compute_models}
\alias{estimate_undetected2molpercell}
\alias{estimate_missingdata}
\alias{compute_genewise_zeroinflation}
\alias{counts2mpc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates noise in single cell data.
}
\description{
Estimates the drop-out model and technical variance from the ERCC spike-ins present in the sample.
}
\usage{
estimate_noiseparameters(spike_counts.df, endogenous_counts.df, spike_conc.df,\cr
  plot = FALSE, sd_inflate = 0, granularity = 300, write.noise.model = TRUE,\cr
  file = "noise_estimation", model_view = c("Observed", "Optimized"),\cr
  total_sampling = 2500, dropout_inflate = 1, alpha_granularity = 0.005)\cr
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spike_counts.df}{
A data frame of spike-in counts for each sample.
}
  \item{endogenous_counts.df}{
A data frame of endogenous gene expression counts for each sample.
}
  \item{spike_conc.df}{
A data frame with spike-in IDs as rownames and expected spike-in counts in the next column.
}
  \item{plot}{
When \code{plot=TRUE} produces plots to investigate quality of data fits with root file name set by \code{file} option.
}
  \item{sd_inflate}{
An optional parameter to modulate the estimated noise. The estimated standard deviation of spike-ins can be scaled by this factor. We recommend leaving the value at 1.
}
  \item{granularity}{
The parameter determines the number of bins for comparison of the quality of fit between the mixed-model and observed data for each alpha. This should be set lower for small datasets and higher for datasets with more observations
}
  \item{write.noise.model}{
When \code{write.noise.model=TRUE} outputs two tab-delimited files containing the dropout effects and noise model parameters; this allows users to apply the noise generation on a seperate high compute node. The root file name is set by \code{file} option.
}
  \item{file}{
Describes the root name for files written out by \code{write.noise.model} and \code{plot} options.
}
  \item{model_view}{
\code{model_view=c("Observed", "Optimized", "Poisson", "Neg. Binomial"} determines the statistical distributions that should be plotted for the ERCC plots output by \code{plot=TRUE}.
}
  \item{total_sampling}{
Defines the number of observations sampled from the noise model in order to compare with actual observations for optimizing alpha.
}
  \item{dropout_inflate}{
A scaling parameter for increasing explicitly the number of drop-outs present beyond those estimated by spike-ins.
}
  \item{alpha_granularity}{
Parameter determines the resolution of alpha values tested for maximum empirical fit to spike-ins.
}
}
\details{
BEARscc consists of three steps: modelling technical variance based on spike-ins (Step 1); simulating technical replicates (Step 2); and clustering simulated replicates (Step 3). In Step 1, an experiment-specific model of technical variability (“noise”) is estimated using observed spike-in read counts. This model consists of two parts. In the first part, expression-dependent variance is approximated by fitting read counts of each spike-in across cells to a mixture model (see Methods). The second part, addresses drop-out effects. Based on the observed drop-out rate for spike-ins of a given concentration, the ‘drop-out injection distribution’ models the likelihood that a given transcript concentration will result in a drop-out. The ‘drop-out recovery distribution’ is estimated from the drop-out injection distribution using Bayes’ theorem and models the likelihood that a transcript that had no observed counts in a cell was a false negative. This funciton performs the first step of BEARscc. For further algoirthmic detail please refer to our manuscript methods.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
\packageAuthor{BEARscc}

Maintainer: \packageMaintainer{BEARscc}
}

\section{Subfunctions of note}{
  \code{estimate_noiseparameters} relies on the following subfunctions to estimate noise. These functions share many common options with the user function. For those options that are internal to the programming; these are annotated to give an idea of flow. For further detail please examine source code in R directory of this package:\cr
  \itemize{
    \item spikes_prepared <- \code{prepare_data(spike_counts.df, spike_conc.df)}
    \item mu2sigma <- \code{estimate_mu2sigma(spikes.prepared, spike_conc.df), plot, sd_inflate, file)}
    \item parameters <- \code{compute_alpha(spikes_prepared, sd_inflate, plot, granularity, file, alpha_granularity)}
    \item   spikes_melted <- melt(data.table(spikes_prepared, keep.rownames = TRUE), id.vars = c("rn", "transcripts"))
    \item models.dt <- \code{compute_models(spikes_melted, parameters)}
    \item counts2mpc.fit <- \code{counts2mpc(spikes_prepared,plot,file)}
    \item \code{estimate_missingdata(spikes_prepared, endogenous_counts.df, counts2mpc.fit, plot, file, dropout_inflate)}
    \item \code{compute_genewise_zeroinflation(endogenous_counts.df)}
    \item \code{spikes.dt <- data.table(spikes_prepared, keep.rownames = TRUE)}
    \item \code{spikes.dt <- spikes.dt[spikes.dt[,rowSums(.SD==0)>1, .SD=c(2:dim(spikes.dt)[2])],][transcripts>0,]}
    \item \code{undetected2mpc <- estimate_undetected2molpercell(spikes.dt, plot, file)}
  }
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(BEARscc_examples)

#For excecution on local machine
estimated_noise <- estimate_noiseparameters(ERCC.counts.df,
                                    data.counts.df,
                                    ERCC.meta.df,
                                    granularity=30,
                                    write.noise.model=FALSE)
estimated_noise

#To save results as files for high power cluster computation later
estimate_noiseparameters(ERCC.counts.df,
                                    data.counts.df,
                                    ERCC.meta.df,
                                    granularity=30,
                                    write.noise.model=TRUE,
                                    file="noise_estimation",
                                    model_view=c("Observed","Optimized"))

}
\keyword{ distribution }% use one of  RShowDoc("KEYWORDS")
\keyword{ models }
