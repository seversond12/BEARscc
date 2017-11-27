\name{counts2mpc}
\alias{counts2mpc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
counts2mpc(x, plot, file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
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
function (x, plot, file) 
{
    x <- melt(data.table(x, keep.rownames = TRUE)[, -1, with = FALSE], 
        id.vars = c("transcripts"))
    x <- x[transcripts > 0.9, ][, mean(value), by = transcripts]
    counts2mpc.fit <- lm(log2(x$V1 + 0.025) ~ log2(x$transcripts))
    if (plot == TRUE) {
        g <- ggplot(data = x, aes(x = log2(transcripts), y = log2(V1 + 
            0.025))) + geom_point() + xlab("\nERCC molecules per cell, log2(molecules)") + 
            ylab("Observed counts, log2(counts) \n") + ggtitle("Correlation of measured counts as a function of actual transcript number.\n") + 
            geom_abline(slope = coef(counts2mpc.fit)[2], intercept = coef(counts2mpc.fit)[1], 
                color = "red")
        ggsave(paste(file, "cor_counts2moleculespercell.pdf", 
            sep = "_"), plot = g, device = "pdf", width = 10.5, 
            height = 7.5, units = "in")
    }
    counts2mpc.fit
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
