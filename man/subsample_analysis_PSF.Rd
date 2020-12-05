\name{subsample_analysis_PSF}
\alias{subsample_analysis_PSF}
\title{Estimation of the probability of detecting new SV using the PSF model.}

\description{Estimation of the probability of detecting new SV using the PSF model.}

\usage{ subsample_analysis_PSF(sv_dataset, subsample) }

\arguments{

\item{sv_dataset}{Table that contains in the 
rows: non overlaping windows in which the genotype-type is described (as a single integer) 
columns: individuals of a population.
It is simlar to vcf file where the columns are the individuals and the rows the loci.}
\item{sv_dataset}{Vector of integer(s) that contain the size of the sample. It is sub-sample as it is required to be smaller thant the input dataset (sv_dataset)  }

}

\examples{
subsample_analysis_PSF(sv_dataset=NULL, subsample=c(20))
}

