% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dna_reads_normalization.R
\name{compute_ploidy}
\alias{compute_ploidy}
\title{Identify ploidy for samples}
\usage{
compute_ploidy(reads, clusters, baseline_cluster)
}
\arguments{
\item{reads}{dna read counts matrix}

\item{clusters}{vector with labels for which cluster each cell belongs to}

\item{baseline_cluster}{one cluster label assumed to be the normal cell population with ploidy 2.}
}
\value{
normalized read counts centered around ploidy 2. amplicons with low read counts are removed
}
\description{
Identify ploidy for samples
}
