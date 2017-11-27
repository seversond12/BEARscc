\name{prepare_probabilities}
\alias{prepare_probabilities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
HPC function analagous to \code{create_noiseinjected_counts}
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
prepare_probabilities(x, probs4detection, parameters, HPC_genewise_permute_count = HPC_genewise_permute_count, HPC_permute_count = HPC_permute_count, HPC_randomizer = HPC_randomizer, total_sampling = 2500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{probs4detection}{
%%     ~~Describe \code{probs4detection} here~~
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
}
  \item{HPC_genewise_permute_count}{
%%     ~~Describe \code{HPC_genewise_permute_count} here~~
}
  \item{HPC_permute_count}{
%%     ~~Describe \code{HPC_permute_count} here~~
}
  \item{HPC_randomizer}{
%%     ~~Describe \code{HPC_randomizer} here~~
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
function (x, probs4detection, parameters, HPC_genewise_permute_count = HPC_genewise_permute_count,
    HPC_permute_count = HPC_permute_count, HPC_randomizer = HPC_randomizer,
    total_sampling = 2500)
{
    probs4detection.genes <- t(data.frame(probs4detection, row.names = "k")[,
        4:eval(dim(probs4detection)[2] - 1)])
    rownames(probs4detection.genes) <- gsub("^X([0-9])", "\\1",
        rownames(probs4detection.genes))
    probs4detection.k <- data.frame(probs4detection[, 2:4, with = FALSE],
        row.names = "k")
    x[, parApply(cl, .SD, 1, HPC_genewise_permute_count, probs4detection.k = probs4detection.k,
        probs4detection.genes = probs4detection.genes, total_sampling = total_sampling,
        HPC_permute_count = HPC_permute_count, HPC_randomizer = HPC_randomizer,
        parameters = parameters)]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
