\name{estimate_missingdata}
\alias{estimate_missingdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimate_missingdata(x, y, counts2mpc.fit, plot, file, dropout_inflate = dropout_inflate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{counts2mpc.fit}{
%%     ~~Describe \code{counts2mpc.fit} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{dropout_inflate}{
%%     ~~Describe \code{dropout_inflate} here~~
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
function (x, y, counts2mpc.fit, plot, file, dropout_inflate = dropout_inflate) 
{
    bayes_list <- compute_genewise_zeroinflation(y)
    ERCC.dt <- data.table(x, keep.rownames = TRUE)
    ERCC.dt <- ERCC.dt[ERCC.dt[, rowSums(.SD == 0) > 1, .SD = c(2:dim(ERCC.dt)[2])], 
        ][transcripts > 0, ]
    undetected2mpc <- estimate_undetected2molpercell(ERCC.dt, 
        plot = plot, file = file)
    kmax <- ceiling(2^(-as.numeric(undetected2mpc[2])/((1/as.numeric(dropout_inflate)) * 
        as.numeric(undetected2mpc[1]))))
    k <- seq(from = 0, to = kmax, by = 1)
    y0.df <- data.frame(k = k)
    y0.df$Py0givenk <- (1/as.numeric(dropout_inflate)) * as.numeric(undetected2mpc[1]) * 
        log2(y0.df$k) + as.numeric(undetected2mpc[2])
    y0.dt <- data.table(y0.df, keep.rownames = TRUE)
    y0.dt[y0.dt$Py0givenk > 1 | y0.dt$Py0givenk == Inf, ]$Py0givenk <- 1
    y0.dt[y0.dt$Py0givenk < 0, ]$Py0givenk <- 0
    y0.dt <- y0.dt[y0.dt$Py0givenk > 0, ]
    y0.dt$Counts <- 2^(coef(counts2mpc.fit)[2] * log2(y0.dt$k) + 
        coef(counts2mpc.fit)[1])
    size <- dim(y0.dt)[1] - 1
    colnames(bayes_list$genewise) <- "py0"
    bayes_list$genewise$pk0 <- (bayes_list$genewise$py0 - sum(y0.dt$Py0givenk[-1])/size)/(1 - 
        sum(y0.dt$Py0givenk[-1])/size)
    if (length(bayes_list$genewise[bayes_list$genewise$pk0 < 
        0, ]$pk0) > 0) {
        bayes_list$genewise[bayes_list$genewise$pk0 < 0, ]$pk0 <- 0
        bayes_list$genewise[bayes_list$genewise$pk0 == 0 & !(bayes_list$genewise$py0 == 
            0), ]$py0 <- sum(y0.dt$Py0givenk[-1])/size
    }
    bayes_list$genewise$pki <- (1 - bayes_list$genewise$pk0) * 
        1/size
    py0givenk <- as.vector(y0.dt$Py0givenk[-1])
    gene.pki <- as.vector(bayes_list$genewise$pki)
    names(gene.pki) <- as.character(rownames(bayes_list$genewise))
    bayes_list$inferred_prob <- py0givenk \%*\% t(gene.pki)
    bayes_list$inferred_prob <- sweep(rbind(t(bayes_list$genewise$pk0), 
        bayes_list$inferred_prob), 2, bayes_list$genewise$py0, 
        "/")
    bayes_list$inferred_prob <- data.table(data.frame(bayes_list$inferred_prob), 
        keep.rownames = TRUE)
    bayes_list$bayes_parameters <- y0.dt
    setkey(bayes_list$bayes_parameters, rn)
    setkey(bayes_list$inferred_prob, rn)
    bayes_list$bayes_parameters <- bayes_list$bayes_parameters[bayes_list$inferred_prob, 
        ]
    setkey(bayes_list$bayes_parameters, k)
    bayes_list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
