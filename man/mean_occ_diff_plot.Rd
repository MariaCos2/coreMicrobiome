% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mean_occ_diff_plot.R
\name{mean_occ_diff_plot}
\alias{mean_occ_diff_plot}
\title{mean_occ_diff_plot}
\usage{
mean_occ_diff_plot(
  otu,
  sample,
  min_abun,
  max_abun,
  abun_step,
  min_thre,
  max_thre,
  thre_step,
  sample_name,
  sample_group
)
}
\arguments{
\item{otu}{a community count data matrix with samples in rows and OTUs/taxa in column,}

\item{sample}{a sample information data.frame. The rownames must match the sample names in the otu,}

\item{max_abun}{the maximum abundance of the loop,}

\item{abun_step}{the abundance step of the loop,}

\item{max_thre}{the maximum threshold of the loop,}

\item{thre_step}{the threshold step of the loop,}

\item{sample_name}{the name of a column that containing the sample name,}

\item{sample_group}{the name of a column that containing the sample type/group,}

\item{mini_abun}{the minimum abundance of the loop,}

\item{mini_thre}{the minimum threshold of the loop,}
}
\value{
A heatmap showing the mean occurrence frequency difference between core and non-core otu under different minimal abundance and threshold.
}
\description{
Plots a heatmap showing the mean occurrence frequency difference between core and non-core otu under different minimal abundance and threshold.
}
\details{
Plots a heatmap showing the mean occurrence frequency difference between core and non-core otu under different minimal abundance and threshold.
}
\examples{
 mean_occ_diff_plot(otu, sample, min_abun = 0, max_abun = 400, abun_step = 100, min_thre = 0, max_thre = 0.1, thre_step = 0.02, sample_name, sample_group)
}
