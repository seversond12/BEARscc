\name{compute_models}
\alias{compute_models}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subfunction in \code{estimate_noiseparameters()}.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute_models(x, parameters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
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
function (x, parameters)
{
    mu2sigma.slope <- parameters$mu2sigma.slope
    sd_inflate <- parameters$sd.inflate
    max.res <- parameters$max.res
    mu2sigma.intercept <- parameters$mu2sigma.intercept
    alpha2mu.slope <- parameters$alpha2mu.slope
    alpha2mu.intercept <- parameters$alpha2mu.intercept
    x2 <- data.frame(rn = "", transcripts = -1, variable = "",
        value = -1)
    for (i in unique(x$rn)) {
        x1 <- x[rn == i]
        t <- x[rn == i]$transcripts
        y <- length(x1$value)
        p <- alpha2mu.slope * log2(mean(x1$value) + 0.025) +
            alpha2mu.intercept
        if (p > 1) {
            p = 1
        }
        if (p < 0) {
            p = 0
        }
        hx.1 <- rnbinom(ceiling(y * p), size = mean(x1$value)^2/((2^(mu2sigma.slope *
            log2(mean(x1$value)) + mu2sigma.intercept + max.res *
            sd_inflate))^2 - mean(x1$value)), mu = mean(x1$value))
        hx.2 <- rpois(ceiling(y * (1 - p)), mean(x1$value))
        hx1 <- rnbinom(y, size = mean(x1$value)^2/((2^(mu2sigma.slope *
            log2(mean(x1$value)) + mu2sigma.intercept + max.res *
            sd_inflate))^2 - mean(x1$value)), mu = mean(x1$value))
        hx2 <- rpois(y, mean(x1$value))
        hx <- sample(append(hx.1, hx.2), y, replace = FALSE)
        sims.m <- data.frame(rn = i, transcripts = t, variable = "Mixed",
            value = hx)
        sims.nb <- data.frame(rn = i, transcripts = t, variable = "NBinom",
            value = hx1)
        sims.p <- data.frame(rn = i, transcripts = t, variable = "Poisson",
            value = hx2)
        x3 <- rbind(x2, sims.m, sims.nb, sims.p, x1)
        x2 <- x3
    }
    x2[-1, ]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
