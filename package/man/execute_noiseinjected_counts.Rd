\name{execute_noiseinjected_counts}
\alias{execute_noiseinjected_counts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subfunction in \code{create_noiseinjected_counts()}.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
execute_noiseinjected_counts(n, noise_parameters, total_sampling)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{noise_parameters}{
%%     ~~Describe \code{noise_parameters} here~~
}
  \item{total_sampling}{
%%     ~~Describe \code{total_sampling} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
%%  ~~who you are~~
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
function (n, noise_parameters, total_sampling)
{
    force(total_sampling)
    print(paste("Creating a noise-injected counts matrix: ",
        n, ".", sep = ""))
    probs4detection.genes <- t(data.frame(noise_parameters$bayes_parameters,
        row.names = "k")[, 4:eval(dim(noise_parameters$bayes_parameters)[2] -
        1)])
    probs4detection.k <- data.frame(noise_parameters$bayes_parameters[,
        2:4, with = FALSE], row.names = "k")
    noisy_counts <- data.table(noise_parameters$original.counts,
        keep.rownames = TRUE)[, apply(.SD, 1, genewise_permute_count,
        probs4detection.k = probs4detection.k, probs4detection.genes = probs4detection.genes,
        parameters = noise_parameters$ERCC_parameters, total_sampling = total_sampling)]
    noisy_counts <- t(noisy_counts)
    rownames(noisy_counts) <- rownames(noise_parameters$original.counts)
    colnames(noisy_counts) <- colnames(noise_parameters$original.counts)
    noisy_counts.dt <- data.table(noisy_counts, keep.rownames = TRUE)
    colnames(noisy_counts.dt)[1] <- "GENE_ID"
    noisy_counts.dt
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
