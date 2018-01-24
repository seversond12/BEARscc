\name{estimate_noiseparameters}
\alias{estimate_noiseparameters}
\alias{prepare_data}
\alias{melt_spikeins}
\alias{plot_mu2sigma}
\alias{estimate_mu2sigma}
\alias{plot_obs2actual}
\alias{iterate_spikeins}
\alias{iterate_alphas}
\alias{plot_alpha2mu}
\alias{compute_alpha}
\alias{subcompute_sample_models}
\alias{sample_models}
\alias{cleanup_model_names}
\alias{subplot_spikein_fits}
\alias{plot_spikein_fits}
\alias{estimate_undetected2molpercell}
\alias{counts2mpc}
\alias{apply_bayes}
\alias{compute_genewise_dropouts}
\alias{build_dropoutmodel}
\alias{create_null_dropout_model}
\alias{estimate_missingdata}
\alias{write_noise_model}
\alias{transcripts}
\alias{V1}

\title{
Estimates noise in single cell data.
}
\description{
Estimates the drop-out model and technical variance
    from spike-ins present in the sample.
}
\usage{
estimate_noiseparameters(SCEList, plot=FALSE, sd_inflate=0, max_cumprob=0.9999,
    bins=10, write.noise.model=TRUE, file="noise_estimation",
    dropout_inflate=1, model_view=c("Observed", "Optimized"),
    alpha_resolution=0.005, tie_function="maximum")
}

\arguments{
\item{SCEList}{
A \code{SingleCellExperiment} object that must contain the observed counts
    matrix as \code{"observed_expression"} in \code{assays}, and must
    have the relevant spike-in samples identified using \code{isSpike()}
    as well as contain the expected actual concentrations of these spike-ins
    as \code{spikeConcentrations} in \code{metadata}. Please see the vignette
    for more detail about constructing the appropriate \code{SCEList}.
}
\item{plot}{
When \code{plot=TRUE} produces plots to investigate quality of
    data fits with root file name set by \code{file} option.
}
\item{sd_inflate}{
An optional parameter to modulate the estimated noise. The estimated
    standard deviation of spike-ins can be scaled by this factor.
    We recommend leaving the value at the default of 0.
}
\item{bins}{
The parameter determines the number of bins for comparison of the quality
    of fit between the mixed-model and observed data for each spike-in
    alpha in order to calculate the relationship between alpha and mean in
    the noise model. This should be set lower for small datasets and higher for
    datasets with more observations
}
\item{max_cumprob}{
Because a cumulative distribution will range from n=0 to a countable
    infinity, the event space needs to be set to cover a reasonable
    fraction of the probability density. This parameter determines the
    the fraction of probability density covered by the event space, which
    in turn defines the highes count number in the event space. We recommend
    users use the default value of 0.9999.
}
\item{write.noise.model}{
When \code{write.noise.model=TRUE} outputs two tab-delimited files
    containing the dropout effects and noise model parameters; this allows
    users to apply the noise generation on a seperate high compute node.
    The root file name is set by \code{file} option.
}
\item{file}{
Describes the root name for files written out by \code{write.noise.model}
    and \code{plot} options.
}
\item{dropout_inflate}{
A scaling parameter for increasing explicitly the number of drop-outs
    present beyond those estimated by spike-ins. The value must be
    greater than 0 or an error will occur. Values below one will
    diminish drop-outs in simulated replicates, and values above
    one will increase drop-outs in simulated replicates. We recommend
    users use the default value of 1.
}
\item{model_view}{
\code{model_view=c("Observed", "Optimized", "Poisson", "Neg. Binomial"}
    determines the statistical distributions that should be plotted for
    the ERCC plots output by \code{plot=TRUE}.
}
\item{alpha_resolution}{
Because the alpha parameter is enumerated discretely and empirically
    evaluated for each value for each spike-in, it is necessary to specify
    the resolution (how small the step is between each explicit alpha test);
    this parameter defines the resolution of alpha values tested for
    maximum empirical fit to spike-ins. It is recommended that users
    utilize the default resolution.
}
\item{tie_function}{
The parameter \code{tie_function=c("minimum", "maximum")} tells
    BEARscc how to handle a tie alpha value for fitting the
    mixture model to an individual spike-in. If \code{maximum},
    then BEARscc will chose the maximum alpha value with the
    best fit; conversely, if \code{minimum} is set, then
    BEARscc will choose the minimum alpha value with the best fit.
}
}

\details{
BEARscc consists of three steps: modelling technical variance based on
    spike-ins (Step 1); simulating technical replicates (Step 2);
    and clustering simulated replicates (Step 3). In Step 1,
    an experiment-specific model of technical variability ("noise")
    is estimated using observed spike-in read counts. This model consists
    of two parts. In the first part, expression-dependent variance is
    approximated by fitting read counts of each spike-in across cells to a
    mixture model (see Methods). The second part, addresses drop-out effects.
    Based on the observed drop-out rate for spike-ins of a given
    concentration, the 'drop-out injection distribution' models the
    likelihood that a given transcript concentration will result in a
    drop-out. The 'drop-out recovery distribution' is estimated from the
    drop-out injection distribution using Bayes' theorem and models the
    likelihood that a transcript that had no observed counts in a cell
    was a false negative. This function performs the first step of BEARscc.
    For further algorithmic detail please refer to our manuscript methods.
}
\value{
The resulting output of \code{estimate_noiseparameters()} is another
\code{SingleCellExperiment} class object; however four new annotations that
describe the drop-out and variance models computed by \code{BEARscc}
have been added to the \code{metadata} of the \code{SingleCellExperiment}
object. Specifically.

\item{\code{dropout_parameters}}{
A \code{data.frame} listing gene-wise parameters necessary for computing
    drop-oout recovery and injection probabilities in order to define the two
    drop-out models for zero observation and positive values within the
    drop-out range by \code{simulate_replicates()}.
}
\item{\code{spikein_parameters}}{
A \code{data.frame} of the estimated noise model parameters utilized by
    \code{simulate_replicates()} to simulate replicates in non-zero
    observations.
}
\item{\code{genewiseDropouts}}{
A \code{data.frame} of the estimated probabilities used in the Bayes'
    calculation of the probabilities described in \code{dropout_parameters}.
    While these are not use in further analysis, they are supplied here for
    the user's reference.
}
}

\author{
\packageAuthor{BEARscc}

Maintainer: \packageMaintainer{BEARscc}
}

\note{
Frequently, the user will want to compute simulated technical replicates in
    a high performance computational environment. While the function outputs
    the necessary information for \code{create_noiseinjected_counts()}, with
    the option \code{write.noise.model=TRUE} users are able to save two tab
    delimited files necessary to run \code{HPC_generate_noise_matrices.R} on
    a high performance computational cluster. The option \code{file} is used
    to indicate the desired root label of the files, "*_bayesianestimates.xls"
    and "*_parameters4randomize.xls".

In the examples section, the parameter, \code{alpha_resolution} is set
    to 0.25, which is a terrible resolution for estimating noise,
    but allows the example to run in reasonable to time for checking
    the help files. We recommend the default parameter:
    \code{alpha_resolution}=0.005.
}

\examples{
library("SingleCellExperiment")
data("BEARscc_examples")

#For execution on local machine
BEAR_examples.sce <- estimate_noiseparameters(BEAR_examples.sce,
    alpha_resolution=0.25, write.noise.model=FALSE)
BEAR_examples.sce

#To save results as files for abnalysis on a
#high performance computational cluster
estimate_noiseparameters(BEAR_examples.sce, write.noise.model=TRUE,
    alpha_resolution=0.25, file="noise_estimation",
    model_view=c("Observed","Optimized"))
}

\keyword{ distribution }
\keyword{ models }
