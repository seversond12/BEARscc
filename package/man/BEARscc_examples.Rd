\name{BEARscc_examples}
\alias{BEARscc_examples}
\alias{data.counts.df}
\alias{ERCC.counts.df}
\alias{ERCC.meta.df}
\docType{data}
\title{
Example data for BEARscc.
}
\description{
A toy dataset for applying BEARscc functions as described in the README on 
    https://bitbucket.org/bsblabludwig/bearscc.git.
}
\usage{data("BEARscc_examples")}

\value{
A long list as described below in the format section.
}

\format{
And R data file with: \cr
\tabular{ccl}{
    df \tab "data.counts.df" \tab Toy dataframe of endogenous counts 
        with 117 genes and 50 samples.\cr
    \cr
    df \tab "ERCC.counts.df" \tab Toy dataframe of ERCC counts with 57 
        spike-ins and 50 samples.\cr
    \cr
    df \tab "ERCC.meta.df" \tab Toy dataframe of spike-in concentration 
        values and spike-in labels as row names..\cr
    \cr
    \tab [,1] \tab Spike-in actual concentration.
}
}

\source{
These data are a subset of observations Drs. Michael White and Richard Owen 
    in the Xin Lu Lab. Samples were sequenced by the Wellcome Trust Center 
    for Genomics, Oxford, UK. These data are available in full with GEO 
    accession number, GSE95155.
}
\references{
Source code and README: <https://bitbucket.org/bsblabludwig/bearscc/overview>
Associated preprint: <https://www.biorxiv.org/content/early/2017/06/05/118919>
}

\examples{
data(BEARscc_examples)
}
\keyword{datasets}
