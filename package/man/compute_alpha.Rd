\name{compute_alpha}
\alias{compute_alpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subfunction in \code{estimate_noiseparameters()}.
}
\description{
Subfunction in \code{estimate_noiseparameters()}.
}
\usage{
compute_alpha(x, estimate_mu2sigma, sd_inflate, plot, granularity, file, alpha_granularity = alpha_granularity)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{estimate_mu2sigma}{
%%     ~~Describe \code{estimate_mu2sigma} here~~
}
  \item{sd_inflate}{
%%     ~~Describe \code{sd_inflate} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{granularity}{
%%     ~~Describe \code{granularity} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{alpha_granularity}{
%%     ~~Describe \code{alpha_granularity} here~~
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
function (x, estimate_mu2sigma, sd_inflate, plot, granularity,
    file, alpha_granularity = alpha_granularity)
{
    mu2sigma <- estimate_mu2sigma(x, plot, sd_inflate, file = file)
    x1 <- melt(data.table(x, keep.rownames = TRUE), id.vars = c("rn",
        "transcripts"))
    if (plot == TRUE) {
        testing <- melt(data.table(x, keep.rownames = TRUE)[,
            -1, with = FALSE], id.vars = c("transcripts"))
        g <- ggplot(data = testing[transcripts >= 1, ], aes(x = factor(round(log10(transcripts +
            0.25), digits = 3)), y = log2(value + 0.25))) + geom_violin(scale = "width") +
            xlab("\nERCC molecules per cell / log10") + ylab("Measured counts / log2\n") +
            ggtitle("The count measurement distribution as a function of actual transcript number.\n")
        ggsave(paste(file, "countsdistribution_vs_moleculespercell.pdf",
            sep = "_"), plot = g, device = "pdf", width = 10.5,
            height = 7.5, units = "in")
    }
    y <- data.table(x1, key = "rn")[, mean(value), by = c("rn",
        "transcripts")]
    setkey(x1, rn)
    colnames(y) <- c("rn", "transcripts", "Mean")
    y <- y[Mean > 0, ]
    ERCCnames <- as.character(y$rn)
    setkey(y, "rn")
    z2 <- data.frame(as.character("DUMMY"), 1, 1, 1, 1)
    colnames(z2) <- c("rn", "Mean", "p.nbinom", "p.poisson",
        "r.value")
    for (j in seq(from = 0, to = 1, by = alpha_granularity)) {
        if (j \%in\% seq(from = 0, to = 1, by = 0.05)) {
            print(paste("Estimating error for ERCCs with alpha = ",
                j, sep = ""))
        }
        k = 1 - j
        y3 <- data.frame(as.character("DUMMY"), 1, 1, 1, 1)
        colnames(y3) <- c("rn", "Mean", "p.nbinom", "p.poisson",
            "r.value")
        for (i in ERCCnames) {
            x2 <- x1[rn == i]
            p <- dim(x1[rn == i])[1]
            y.name <- as.character(y[rn == i]$rn)
            y.mean <- y[rn == i]$Mean
            hx1 <- rnbinom(p * j, size = mean(x2$value)^2/((2^(mu2sigma[1,
                1] * log2(mean(x2$value)) + mu2sigma[1, 2] +
                mu2sigma[1, 3]))^2 - mean(x2$value)), mu = mean(x2$value))
            hx2 <- rpois(p * k, mean(x2$value))
            if (j == 0) {
                hx.mixed <- hx2
            }
            else if (j == 1) {
                hx.mixed <- hx1
            }
            else {
                hx.mixed <- append(hx1, hx2)
            }
            brakes <- seq(from = 0, to = max(hx.mixed, x2$value) +
                max(hx.mixed, x2$value)/granularity, by = max(hx.mixed,
                x2$value)/granularity)
            h1 <- hist(x2$value, plot = FALSE, breaks = brakes)
            h.mixed <- hist(hx.mixed, plot = FALSE, breaks = brakes)
            rvalue <- sum(abs(h.mixed$counts - h1$counts)/(2 *
                p))
            y2 <- data.frame(rn = y.name, Mean = y.mean, p.nbinom = j,
                p.poisson = k, r.value = rvalue)
            y4 <- rbind(y3, y2)
            y3 <- y4
        }
        z1 <- rbind(z2, y3)
        z2 <- data.table(z1, key = "rn")[!(rn == "DUMMY"), ]
    }
    z2 <- data.table(z2, key = "rn")
    DT <- z2[, min(r.value), by = rn]
    y1 <- data.frame(as.character("DUMMY"), 1, 1, 1, 1)
    colnames(y1) <- c("rn", "Mean", "p.nbinom", "p.poisson",
        "r.value")
    for (i in DT$rn) {
        y2 <- z2[rn == i & r.value == DT[rn == i]$V1, ]
        y3 <- rbind(y1, y2)
        y1 <- y3
    }
    rvalues.dt <- data.table(y1[-1], key = "rn")
    rvalues.dt <- rvalues.dt[rvalues.dt$Mean >= 1, ]
    alpha2mean <- lm(rvalues.dt$p.nbinom ~ log2(rvalues.dt$Mean +
        0.025))
    if (plot == TRUE) {
        g <- ggplot(data = rvalues.dt, aes(x = log2(Mean + 0.025),
            y = p.nbinom)) + geom_point(fill = "red", color = "black",
            size = 2.2, alpha = 0.6, pch = 21) + xlab("\nObserved mean ERCC expression, log2(expression)") +
            ylab("Neg. binomial contribution parameter, alpha\n") +
            ggtitle("The neg. binomial contribution is a function of gene expression.\n") +
            geom_abline(slope = coef(alpha2mean)[2], intercept = coef(alpha2mean)[1])
        ggsave(paste(file, "alpha2mu_correlation.pdf", sep = "_"),
            plot = g, device = "pdf", width = 10.5, height = 7.5,
            units = "in")
    }
    alpha2mu <- data.frame(alpha2mu.slope = coef(alpha2mean)[2],
        alpha2mu.intercept = coef(alpha2mean)[1], sd.inflate = sd_inflate)
    parameters <- data.frame(cbind(mu2sigma, alpha2mu), row.names = c("parameter.value"))
    parameters
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
