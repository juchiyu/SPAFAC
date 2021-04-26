% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sDICA.R
\name{sparseDiCA}
\alias{sparseDiCA}
\title{Group sparse Discriminant Correspondence Analysis}
\usage{
sparseDiCA(
  DATA,
  design,
  components = 0,
  tol = .Machine$double.eps,
  make_data_nominal = TRUE,
  make_design_nominal = TRUE,
  doublecentering = TRUE,
  init = "svd",
  initLeft = NULL,
  initRight = NULL,
  seed = NULL,
  rdsLeft = rep(1, components),
  rdsRight = rep(1, components),
  grpLeft = NULL,
  grpRight = NULL,
  orthogonality = "loadings",
  OrthSpaceLeft = NULL,
  OrthSpaceRight = NULL,
  projPriority = "orth",
  projPriorityLeft = projPriority,
  projPriorityRight = projPriority,
  itermaxALS = 1000,
  itermaxPOCS = 1000,
  epsALS = 1e-10,
  epsPOCS = 1e-10
)
}
\arguments{
\item{DATA}{the contingency table}

\item{design}{the design vector of the observations (i.e., the rows) of DATA}

\item{grpLeft}{the grouping vector for the rows}

\item{grpRight}{the grouping vector for the columns}

\item{epsPOCS}{}

\item{k}{the number of dimensions}
}
\value{

}
\description{
Group sparse Discriminant Correspondence Analysis
}