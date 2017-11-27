\name{estimate_undetected2molpercell}
\alias{estimate_undetected2molpercell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subfunction in \code{estimate_noiseparameters()}.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimate_undetected2molpercell(x, plot, file)
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
    undetected <- data.table(x[, rowSums(.SD == 0), .SD = c(2:(dim(x)[2] -
        1))]/(dim(x)[2] - 2))
    undetected$rn <- as.character(x$rn)
    undetected$transcripts <- as.numeric(x$transcripts)
    undetected <- undetected[undetected$transcripts > 0.9, ]
    undetected2mpc.fit <- lm(undetected$V1 ~ log2(undetected$transcripts))
    undetected2molpercell <- data.frame(undetected2mpc.slope = coef(undetected2mpc.fit)[2],
        undetected2mpc.intercept = coef(undetected2mpc.fit)[1])
    if (plot == TRUE) {
        g <- ggplot(data = undetected, aes(x = log2(transcripts),
            y = V1)) + geom_point() + xlab("\nERCC actual counts, log2(molecules per cell)") +
            ylab("Fraction of cells where observed ERCC count is zero") +
            ggtitle("Conditional probability of observing a zero given transcript concentration\n") +
            geom_abline(slope = coef(undetected2mpc.fit)[2],
                color = "red", intercept = coef(undetected2mpc.fit)[1])
        ggsave(paste(file, "undetected2molpercell.pdf", sep = "_"),
            plot = g, device = "pdf", width = 11, height = 8.5,
            units = "in")
    }
    undetected2molpercell
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
