\name{HPC_randomizer}
\alias{HPC_randomizer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subfunction in \code{prepare_probabilities()}.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HPC_randomizer(x, parameters, total_sampling)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x, parameters, total_sampling)
{
    NB_intercept <- as.numeric(parameters$alpha2mu.intercept)
    NB_slope <- as.numeric(parameters$alpha2mu.slope)
    Mean2SDcor_residual <- as.numeric(parameters$max.res)
    Mean2SDcor_slope <- as.numeric(parameters$mu2sigma.slope)
    Mean2SDcor_intercept <- as.numeric(parameters$mu2sigma.intercept)
    sd_inflate <- as.numeric(parameters$sd.inflate)
    p <- NB_intercept + NB_slope * log2(x + 0.025)
    if (p < 0) {
        p = 0
    }
    if (p > 1) {
        p = 1
    }
    hx1 <- rnbinom(total_sampling * p, size = x^2/((2^(Mean2SDcor_slope *
        log2(x) + Mean2SDcor_intercept + Mean2SDcor_residual *
        sd_inflate))^2 - x), mu = x)
    hx2 <- rpois(total_sampling * (1 - p), x)
    sample(append(hx1, hx2), size = 1, replace = TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
