% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FeatDyn_xCT.R
\name{FeatDyn_cellType}
\alias{FeatDyn_cellType}
\title{Feature dynamics per cell type}
\usage{
FeatDyn_cellType(data, Feature, tp_explicit = TRUE)
}
\arguments{
\item{data}{Entire data.frame}

\item{tp_explicit}{Take into account time_points?}

\item{Features}{The selected feature}
}
\value{
A plot representing the dynamic change of the feature
}
\description{
This function generates a plot representing the dynamic change of the feature per cell type
}
