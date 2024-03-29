% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mppGE_QTL_R2.R
\name{mppGE_QTL_R2}
\alias{mppGE_QTL_R2}
\title{MPP GxE QTL R2}
\usage{
mppGE_QTL_R2(
  mppData,
  trait,
  Q.eff = "cr",
  VCOV = "ID",
  QTL = NULL,
  glb.only = FALSE,
  workspace = 8e+06
)
}
\arguments{
\item{mppData}{An object of class \code{mppData}.}

\item{trait}{\code{Character vector} specifying which traits (environments) should be used.}

\item{Q.eff}{\code{Character} expression indicating the assumption concerning
the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
for ancestral; 4) "biall" for a bi-allelic. Default = "cr".}

\item{VCOV}{VCOV \code{Character} expression defining the type of variance
covariance structure used. "ID" for identity, "CSRT" for within environment
cross-specific residual terms, "CS_CSRT" for compound symmetry with within
environment cross-specific residual terms. Default = "CS_CSRT".}

\item{QTL}{Object of class \code{QTLlist} representing a list of
selected marker positions obtained with the function QTL_select() or
a vector of \code{character} marker positions names. Default = NULL.}

\item{glb.only}{\code{Logical} value. If glb.only = TRUE, only the global and
global adjusted R squared will be returned. Default = FALSE.}

\item{workspace}{Size of workspace for the REML routines measured in double
precision words (groups of 8 bytes). The default is workspace = 8e6.}
}
\value{
Return:

\item{}{}
}
\description{
Compute MPP GxE QTL R2.
}
\examples{

library(asreml)

data(mppData_GE)

Qpos <- c("PZE.105068880", "PZE.106098900")

R2 <- mppGE_QTL_R2(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
                   Q.eff = 'par', QTL = Qpos)

}
\author{
Vincent Garin
}
