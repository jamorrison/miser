% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addAnnotations.R
\name{addAnnotations}
\alias{addAnnotations}
\title{Add annotations to object returned from sesamizeGEO}
\usage{
addAnnotations(x, pheno)
}
\arguments{
\item{x}{Object returned from sesamizeGEO}

\item{pheno}{phenoData data frame from GEOquery}
}
\value{
Object returned from sesamizeGEO with added annotations
}
\description{
GEOmetadb requires a large memory overhead. This function uses the phenoData
data frame to supply the annotation information, rather than using GEOmetadb
like addCharacteristics
}
\examples{

}
