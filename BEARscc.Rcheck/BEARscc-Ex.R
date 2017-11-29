pkgname <- "BEARscc"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "BEARscc-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('BEARscc')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("BEARscc_examples")
### * BEARscc_examples

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: BEARscc_examples
### Title: Example data for BEARscc.
### Aliases: BEARscc_examples data.counts.df ERCC.counts.df ERCC.meta.df
### Keywords: datasets

### ** Examples

data(BEARscc_examples)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("BEARscc_examples", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("analysis_examples")
### * analysis_examples

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: analysis_examples
### Title: BEARscc downstream example objects.
### Aliases: analysis_examples noise_consensus cluster.list sim_replicates
###   estimated_noise recluster clusters.df BEARscc_clusts.df .Random.seed
### Keywords: datasets

### ** Examples

data(analysis_examples)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("analysis_examples", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cluster_consensus")
### * cluster_consensus

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cluster_consensus
### Title: Cluster the consensus matrix.
### Aliases: cluster_consensus
### Keywords: cluster optimize

### ** Examples

data(analysis_examples)

vector <- seq(from=2, to=5, by=1)
BEARscc_clusts.df <- cluster_consensus(noise_consensus, vector)
BEARscc_clusts.df




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cluster_consensus", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_consensus")
### * compute_consensus

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_consensus
### Title: Compute consensus matrix.
### Aliases: compute_consensus create_cm
### Keywords: models error

### ** Examples

data("analysis_examples")

noise_consensus <- compute_consensus(clusters.df)
noise_consensus



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_consensus", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("create_noiseinjected_counts")
### * create_noiseinjected_counts

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: create_noiseinjected_counts
### Title: Computes BEARscc simulated technical replicates.
### Aliases: create_noiseinjected_counts HPC_genewise_permute_count
###   HPC_permute_count HPC_randomizer prepare_probabilities randomizer
###   permute_count genewise_permute_count execute_noiseinjected_counts
### Keywords: models robust

### ** Examples

data(analysis_examples)

sim_replicates<-create_noiseinjected_counts(estimated_noise, n=10)
sim_replicates




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("create_noiseinjected_counts", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("estimate_noiseparameters")
### * estimate_noiseparameters

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: estimate_noiseparameters
### Title: Estimates noise in single cell data.
### Aliases: estimate_noiseparameters prepare_data estimate_mu2sigma
###   compute_alpha compute_models estimate_undetected2molpercell
###   estimate_missingdata compute_genewise_zeroinflation counts2mpc
### Keywords: distribution models

### ** Examples

data(BEARscc_examples)

#For excecution on local machine
estimated_noise <- estimate_noiseparameters(ERCC.counts.df,
                                    data.counts.df,
                                    ERCC.meta.df,
                                    granularity=30,
                                    write.noise.model=FALSE)
estimated_noise

#To save results as files for high performance computation cluster later
estimate_noiseparameters(ERCC.counts.df,
                                    data.counts.df,
                                    ERCC.meta.df,
                                    granularity=30,
                                    write.noise.model=TRUE,
                                    file="noise_estimation",
                                    model_view=c("Observed","Optimized"))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("estimate_noiseparameters", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("report_cell_metrics")
### * report_cell_metrics

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: report_cell_metrics
### Title: Reports BEARscc metrics for cells.
### Aliases: report_cell_metrics calculate_cell_metrics
###   calculate_cell_metrics_by_cluster
### Keywords: list

### ** Examples

data(analysis_examples)

cell_scores.dt <- report_cell_metrics(BEARscc_clusts.df, noise_consensus)
cell_scores.dt



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("report_cell_metrics", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("report_cluster_metrics")
### * report_cluster_metrics

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: report_cluster_metrics
### Title: Reports BEARscc metrics for clusters.
### Aliases: report_cluster_metrics calculate_cluster_metrics
###   calculate_cluster_metrics_by_cluster
### Keywords: list

### ** Examples

data(analysis_examples)

cluster_scores.dt <- report_cluster_metrics(BEARscc_clusts.df,noise_consensus, plot=TRUE, file="example")
clusters.dt



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("report_cluster_metrics", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
