% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/net_viz.R
\name{net_viz}
\alias{net_viz}
\title{Generate, Save, and Plot Networks for Each RT Window}
\usage{
net_viz(
  data,
  image_format = ".jpeg",
  output_xlsx = "network_data.xlsx",
  output_dir = "network_plots"
)
}
\arguments{
\item{data}{A data frame containing network data with at least the columns: \code{Window}, \code{mz1}, \code{mz2}, and \code{molecule_match}.}

\item{image_format}{Specify network image file format from: \code{.jpeg}, \code{.png}, \code{.tiff}, or \code{.svg}.}

\item{output_xlsx}{Character. The name of the output Excel file where network data will be saved. Default is "network_data.xlsx".}

\item{output_dir}{Character. The directory where network plots will be saved. Default is "network_plots".}
}
\value{
No return value. Outputs an Excel file with network data and saves network plots as JPEG images.
}
\description{
This function generates network layouts for each RT window, saves the network data to an Excel file,
and produces network plots that are saved as JPEG files.
}
\examples{
library(igraph)
library(ggraph)
sample_data <- data.frame(
  Window = rep(1:3, each = 5),
  mz1 = runif(15, 100, 500),
  mz2 = runif(15, 100, 500),
  molecule_match = sample(c("A", "B", "C"), 15, replace = TRUE)
)
net_viz(sample_data)

}
