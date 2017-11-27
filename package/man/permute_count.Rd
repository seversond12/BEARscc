\name{permute_count}
\alias{permute_count}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subfunction in \code{create_noiseinjected_counts()}.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
permute_count(x, probs4detection.k, probabilityA, parameters, total_sampling)
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
function (x, probs4detection.k, probabilityA, parameters, total_sampling)
{
    force(parameters)
    if (x == 0) {
        nx <- sample(as.numeric(rownames(probs4detection.k)),
            1, prob = probabilityA, replace = TRUE)
        probabilityB <- probs4detection.k[as.character(ceiling(nx)),
            ]$Py0givenk
        nx <- sample(c(0, nx), 1, prob = as.numeric(c(probabilityB,
            1 - probabilityB)))
        if (nx > 0) {
            nx <- randomizer(probs4detection.k[as.character(nx),
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
                nx <- randomizer(nx, parameters, total_sampling)
            }
        }
        else {
            nx <- randomizer(x, parameters, total_sampling)
        }
    }
    nx
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
