% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_new_data.R
\name{gen_new_data}
\alias{gen_new_data}
\title{Generate new data based on clusters}
\usage{
gen_new_data(data, cluster = "clusters", method = c("mean", "sum"))
}
\arguments{
\item{data}{Your data}

\item{cluster}{Name of the cluster column}

\item{method}{Method applied to generate the new data based on clusters (mean or sum)}
}
\value{
It returns a new dataframe generated based on the clusters
}
\description{
This function generates a new dataframe based on clusters generated by the function "divHier()" (or in the same format), 
it can be generated by the sum or the mean of the rows that are in the same cluster.
}
