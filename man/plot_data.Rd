\name{plot_data}
\alias{plot_data}
\docType{data}
\title{
  Example plot data
}
\description{

  Example plot data from the maize EU-NAM Flint population (Bauer et al. 2013,
  Lehermeier et al. 2014).


}
\usage{data(plot_data)}
\format{
  \code{data.frame}
}
\details{

Sample plot data from the maize EU-NAM Flint population for two environments: La Coruna (CIAM) and Roggenstein (TUM). The 'genotype' column represents the genotype identifiers. The 'check' column is composed of genotype identifier for the genotypes present in the \code{\link{mppData_GE}} object. These genotypes are the one that have marker information that will be used for the QTL detection. These genotypes come from 5 selected crosses. The rest of the check column is composed of genotypes identifiers that are replicated in at least one environment. These genotypes will be used as fixed checks in the one stage analysis QTL detection. The 'cross' column contains the cross indicators. The 'env' column contains the environment indicators. The last columns contain the phenotypic values values. The raw plot data were obtained here: http://www.genetics.org/content/198/1/3/suppl/DC1

The \code{plot_data} argument is generally composed with the following columns: 'genotype', 'check', 'cross', 'env', experimental design variables, traits.

The 'genotype' column contains the genotype (genotype + checks) identifiers.

The 'check' column indicate if the genotype has marker data ('genotype') or if it is a check (same genotype identifier as genotype or 'check').

The 'cross' column indicate to which cross the genotypes belongs. For the check, they can either come from a particular cross or from an 'Extra' cross.

The 'env' column specifies the environment in which the genotype was characterized.

The next columns are used for experimental design variables. The column names must be the same as the one that will be used in 'exp_des_form' in the functions. The traits are numerical values with the same column names as the one used for the 'trait' argument.


}

\references{

  Bauer, E., Falque, M., Walter, H., Bauland, C., Camisan, C., Campo, L., ...
  & Altmann, T. (2013). Intraspecific variation of recombination rate in maize.
  Genome biology, 14(9), R103.

  Lehermeier, C., Krämer, N., Bauer, E., Bauland, C., Camisan, C., Campo, L.,
  ... & Moreau, L. (2014). Usefulness of multiparental populations of maize
  (Zea mays L.) for genome-based prediction. Genetics, 198(1), 3-16.

}

\examples{
  data(plot_data)
}

\seealso{
  \code{\link{mppData_GE}}
}


\keyword{datasets}
