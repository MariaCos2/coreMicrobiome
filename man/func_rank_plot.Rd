% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_rank_plot.R
\name{func_rank_plot}
\alias{func_rank_plot}
\title{func_rank_plot}
\usage{
func_rank_plot(functional_profile, sample, mini_tpm, threshold)
}
\arguments{
\item{functional_profile}{a community count data matrix with samples in rows and OTUs/taxa in column,}

\item{sample}{a sample information data.frame. The rownames must match the sample names in the otu,}

\item{mini_tpm}{a value indicating whether the function is present,}

\item{threshold}{a value indicating the final percent increase in beta-diversity,}
}
\value{
A scatter plot showing the relationship between TPM and ranked otu with color denoting and core or non-core functional otu and a line indicating the threshold and a scatter plot showing the relationship between average occurancy and average TPM with color denoting and core or non-core functional otu and a line indicating the threshold.
}
\description{
Plots TPM vs ranked otu and avg occ vs avg TPM.
}
\details{
Plots TPM vs ranked otu and avg occ vs avg TPM.
}
\examples{
 func_rank_plot(functional_profile, sample, mini_tpm=10000, threshold=0.5)
}
