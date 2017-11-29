\name{create_noiseinjected_counts}
\alias{create_noiseinjected_counts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes BEARscc simulated technical replicates.
}
\description{
Computes BEARscc simulated technical replicates from the previously estimated noise parameters computed with the function \code{estimate_noise_paramters()}.
}
\usage{
create_noiseinjected_counts(noise_parameters, total_sampling = 2500, n = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{noise_parameters}{
%%     ~~Describe \code{noise_parameters} here~~
}
  \item{total_sampling}{
%%     ~~Describe \code{total_sampling} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (noise_parameters, total_sampling = 2500, n = 3)
{
    force(total_sampling)
    vector <- as.character(seq(from = 1, to = n, by = 1))
    vector.names <- gsub("^", "Iteration_", vector)
    vector <- data.table(t(vector))
    colnames(vector) <- vector.names
    noisy_counts.list <- lapply(vector, execute_noiseinjected_counts,
        noise_parameters = noise_parameters, total_sampling = total_sampling)
    noisy_counts.list$Original_counts <- data.table(noise_parameters$original.counts,
        keep.rownames = TRUE)
    colnames(noisy_counts.list$Original_counts)[1] <- "GENE_ID"
    noisy_counts.list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
