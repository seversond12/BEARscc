\name{genewise_permute_count}
\alias{genewise_permute_count}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subfunction in \code{create_noiseinjected_counts()}.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
genewise_permute_count(x, probs4detection.k, probs4detection.genes, parameters, total_sampling)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{probs4detection.k}{
%%     ~~Describe \code{probs4detection.k} here~~
}
  \item{probs4detection.genes}{
%%     ~~Describe \code{probs4detection.genes} here~~
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
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
function (x, probs4detection.k, probs4detection.genes, parameters,
    total_sampling)
{
    probabilityA <- probs4detection.genes[gsub("-", ".", x[1]),
        ]
    force(total_sampling)
    apply(data.frame(as.numeric(x[-1])), 1, permute_count, probs4detection.k,
        probabilityA = probabilityA, parameters, total_sampling)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
