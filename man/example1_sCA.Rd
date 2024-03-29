% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example1_sCA}
\alias{example1_sCA}
\title{USA Deaths in 2018 by Age and Cause
Causes of death in the US in 2019}
\format{
A data frame with 399 rows and 3 columns:
\describe{
\item{AgeGroup}{Age group of death, grouped in 5-year intervals, with a "100+" category.}
\item{CauseOfDeath}{Cause of death, categorized according to the ICD-10.}
\item{NumberOfDeaths}{The number of deaths for each cause at each age range.}
}
}
\source{
International Statistical Classification of Diseases and Related Health Problems (ICD-10).
}
\usage{
example1_sCA
}
\description{
This dataset contains information on the number of deaths in the USA in 2018,
categorized by age and cause of death. The causes of death are classified
according to the International Statistical Classification of Diseases and
Related Health Problems (ICD-10). Ages of death are grouped in 5-year intervals,
starting from 1 year old, with all deaths above 100 years grouped into a "100+"
category. Deaths before the age of 1 are excluded due to their predominance
in the perinatal causes category, which could overshadow other effects.
This exclusion is based on the observation that perinatal causes were rare in other age groups,
with 57 out of 118 cases occurring before 1 year old and the remaining 61 cases
spread across the other 21 age groups. The dataset represents a 21x19
contingency table of counts, showing the number of deaths for each cause
at each age range.
}
\references{
\url{https://www.who.int/classifications/icd/icdonlineversions/en/}
}
\keyword{datasets}
