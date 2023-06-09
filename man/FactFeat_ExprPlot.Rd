% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FactFeat_ExprPlot.R
\name{FactFeat_ExprPlot}
\alias{FactFeat_ExprPlot}
\title{Factor - Feature expression comparison}
\usage{
FactFeat_ExprPlot(Feat_data, Fact_data, Features)
}
\arguments{
\item{Feat_data}{Data.frame of feature expressions}

\item{Fact_data}{Data.frame of factor expression}

\item{Features}{Vector of feature names}
}
\value{
A list containing the plots
}
\description{
This function generates several plots comparing the expression of each feature
and the factor taking into account the parameters time-point and cell-type.
}
