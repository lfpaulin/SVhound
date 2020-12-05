\name{pitman.probability_new_next}
\alias{pitman.probability_new_next}
\title{Compute the probability of observing a new element in the next observation (see PSF).}

\description{Compute the probability of observing a new element in the next observation (see PSF).}

\usage{ pitman.probability_new_next(theta, sigma, k, n) }

\arguments{

\item{theta}{PSF model parameter .} 
\item{sigma}{PSF model parameter .} 
\item{k}{Number of distinct observed elemtens .} 
\item{n}{Total number of observed elements.}

}

\examples{
x <- ewens.probability_new_next(theta=20, sigma=0.001, k=100, n=100000)
}

