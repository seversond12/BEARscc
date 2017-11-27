\name{estimate_mu2sigma}
\alias{estimate_mu2sigma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subfunction in \code{estimate_noiseparameters()}.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimate_mu2sigma(x, plot, sd_inflate, file = "Rplot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{sd_inflate}{
%%     ~~Describe \code{sd_inflate} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
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
function (x, plot, sd_inflate, file = "Rplot")
{
    x1 <- melt(data.table(x, keep.rownames = TRUE), id.vars = c("rn",
        "transcripts"))
    ERCC.sd <- data.table(x1[, sd(value), by = c("rn", "transcripts")])
    ERCC.sd$mean <- x1[, mean(value), by = c("rn", "transcripts")]$V1
    ERCC.sd[ERCC.sd$V1 == 0]$V1 <- NA
    ERCC.sd[ERCC.sd$V1 == 0]$mean <- NA
    mu2sigma.fit <- lm(log2(ERCC.sd$V1) ~ log2(ERCC.sd$mean))
    mu2sigma <- data.frame(mu2sigma.slope = coef(mu2sigma.fit)[2],
        mu2sigma.intercept = coef(mu2sigma.fit)[1], max.res = max(mu2sigma.fit$residuals))
    if (plot == TRUE) {
        g <- ggplot(data = ERCC.sd, aes(x = log2(mean), y = log2(V1))) +
            geom_point() + geom_abline(slope = coef(mu2sigma.fit)[2],
            intercept = coef(mu2sigma.fit)[1]) + xlab("\nMean ERCC expression / log2") +
            ylab("Standard deviation of ERCC expression / log2\n") +
            ggtitle("Relationship between standard deviation and mean of ERCC spike in measurements\n") +
            geom_abline(slope = coef(mu2sigma.fit)[2], color = "red",
                intercept = coef(mu2sigma.fit)[1] + sd_inflate *
                  max(mu2sigma.fit$residuals))
        ggsave(paste(file, "mu2sigma.pdf", sep = "_"), plot = g,
            device = "pdf", width = 10.5, height = 7.5, units = "in")
    }
    mu2sigma
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
