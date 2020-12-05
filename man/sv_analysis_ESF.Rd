\name{sv_analysis_ESF}
\alias{sv_analysis_ESF}
\title{Estimation of the probability of detecting new SV using the ESF model.}

\description{Estimation of the probability of detecting new SV using the ESF model.}

\usage{ sv_analysis_ESF(sv_dataset) }

\arguments{

\item{sv_dataset}{Table that contains in the 
rows: non overlaping windows in which the genotype-type is described (as a single integer) 
columns: individuals of a population.
It is simlar to vcf file where the columns are the individuals and the rows the loci.}

}

\examples{
sv_analysis_ESF(sv_dataset)
}

