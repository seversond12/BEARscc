\name{estimate_noiseparameters}
\alias{estimate_noiseparameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates noise in single cell data.
}
\description{
Estimates the drop-out model and technical variance from the ERCC spike-ins present in the sample.
}
\usage{
estimate_noiseparameters(x, y, z, plot = FALSE, sd_inflate = 0, granularity = 300, write.noise.model = TRUE, file = "noise_estimation", model_view = c("Observed", "Optimized"), total_sampling = 2500, dropout_inflate = 1, alpha_granularity = 0.005)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data frame of spike-in counts for each sample.
}
  \item{y}{
A data frame of endogenous gene expression counts for each sample.
}
  \item{z}{
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, z, plot = FALSE, sd_inflate = 0, granularity = 300,
    write.noise.model = TRUE, file = "noise_estimation", model_view = c("Observed",
        "Optimized"), total_sampling = 2500, dropout_inflate = 1,
    alpha_granularity = 0.005)
{
    force(granularity)
    force(sd_inflate)
    force(alpha_granularity)
    ERCC.prepared.counts <- prepare_data(x, y)
    actual.prepared.counts <- z
    print("Fitting parameter alpha to establish ERCC-derived noise model.")
    parameters <- compute_alpha(ERCC.prepared.counts, estimate_mu2sigma,
        granularity = granularity, sd_inflate = sd_inflate, plot = plot,
        file = file, alpha_granularity = alpha_granularity)
    ERCC.m.counts <- melt(data.table(ERCC.prepared.counts, keep.rownames = TRUE),
        id.vars = c("rn", "transcripts"))
    ERCC.m.valid <- ERCC.m.counts[, mean(value) > 0, by = "rn"]
    setkey(ERCC.m.valid, "rn")
    setkey(ERCC.m.counts, "rn")
    ERCC.m.counts <- ERCC.m.valid[ERCC.m.counts, ][V1 == TRUE,
        ][, -2, with = FALSE]
    models.dt <- compute_models(ERCC.m.counts, parameters)
    models.dt$distribution <- gsub(TRUE, "Observed", models.dt$variable !=
        "Mixed" & models.dt$variable != "NBinom" & models.dt$variable !=
        "Poisson")
    models.dt[models.dt$variable == "Mixed", ]$distribution <- "Optimized"
    models.dt[models.dt$variable == "Poisson"]$distribution <- "Poisson"
    models.dt[models.dt$variable == "NBinom"]$distribution <- "Neg. Binomial"
    if (plot == TRUE) {
        ERCCnames <- as.character(unique(models.dt[order(transcripts)]$rn))
        round(length(ERCCnames)/6)
        dir.create("ERCCfit_plots")
        for (i in seq(from = 1, to = round(length(ERCCnames)/6),
            by = 1)) {
            p <- ERCCnames[seq(from = i, to = length(ERCCnames),
                by = round(length(ERCCnames)/6))]
            g <- ggplot(data = models.dt[rn \%in\% p, ][distribution \%in\%
                model_view], aes(x = value)) + geom_histogram(aes(fill = factor(distribution,
                levels = c("Poisson", "Neg. Binomial", "Optimized",
                  "Observed"))), alpha = 0.3, position = "identity",
                bins = granularity) + facet_wrap(~rn, ncol = 3,
                nrow = 4, scales = "free") + ylab("Observation density\n") +
                labs(fill = "Model") + xlab("\nERCC expression / normalized counts") +
                ggtitle("Approximating technical noise using a refined mixed distribution.\n") +
                scale_color_discrete(guide = FALSE) + theme(axis.text.x = element_text(angle = 90,
                hjust = 1, size = 10))
            savefile <- paste("ERCCfit_plots/", paste(file, "ERCCfits_subset",
                i, "plot.pdf", sep = "_"), sep = "")
            ggsave(savefile, plot = g, device = "pdf", width = 10.5,
                height = 7.5, units = "in")
        }
    }
    counts2mpc.fit <- counts2mpc(ERCC.prepared.counts, plot = plot,
        file = file)
    bayes_dropouts <- estimate_missingdata(ERCC.prepared.counts,
        actual.prepared.counts, counts2mpc.fit, plot = plot,
        file = file, dropout_inflate = dropout_inflate)
    bayes_dropouts$ERCC_parameters <- parameters
    noise_model <- bayes_dropouts
    noise_model$models.dt <- models.dt
    noise_model$original.counts <- z
    noise_model$ERCC.concentrations <- x
    noise_model$ERCC.counts <- y
    if (write.noise.model == TRUE) {
        write.table(data.table(noise_model$ERCC_parameters, keep.rownames = TRUE),
            file = paste(file, "parameters4randomize.xls", sep = "_"),
            quote = FALSE, sep = "\t", row.names = FALSE)
        write.table(noise_model$bayes_parameters, file = paste(file,
            "bayesianestimates.xls", sep = "_"), quote = FALSE,
            sep = "\t", row.names = FALSE)
        print("Parameters have been saved into the current directory in bayesianestimates.xls and parameters4randomize.xls.")
        noise_model
    }
    else {
        noise_model
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
