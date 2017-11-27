\name{HPC_permute_count}
\alias{HPC_permute_count}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subfunction in \code{prepare_probabilities()}.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HPC_permute_count(x, probs4detection.k, probabilityA, parameters, gene, HPC_randomizer, total_sampling)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{probs4detection.k}{
%%     ~~Describe \code{probs4detection.k} here~~
}
  \item{probabilityA}{
%%     ~~Describe \code{probabilityA} here~~
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
}
  \item{gene}{
%%     ~~Describe \code{gene} here~~
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
function (x, probs4detection.k, probabilityA, parameters, gene,
    HPC_randomizer, total_sampling)
{
    force(parameters)
    if (x == 0) {
        nx <- sample(as.numeric(rownames(probs4detection.k)),
            1, prob = probabilityA, replace = TRUE)
        if (nx > 0) {
            nx <- HPC_randomizer(probs4detection.k[as.character(nx),
                ]$Counts, parameters, total_sampling)
        }
    }
    else {
        if (ceiling(x) < max(as.numeric(rownames(probs4detection.k)))) {
            probabilityB <- as.numeric(probs4detection.k[as.character(ceiling(x)),
                ]$Py0givenk)
            nx <- sample(c(0, x), 1, prob = as.numeric(c(probabilityB,
                1 - probabilityB)))
            if (nx > 0) {
                nx <- HPC_randomizer(nx, parameters, total_sampling)
            }
        }
        else {
            nx <- HPC_randomizer(x, parameters, total_sampling)
        }
    }
    nx
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
