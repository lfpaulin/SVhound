\name{ewens.probability_new_next}
\alias{ewens.probability_new_next}
\title{Compute the probability of observing a new element in the next observation (see ESF).}

\description{Compute the probability of observing a new element in the next observation (see ESF).}

\usage{ ewens.probability_new_next(theta, n) }

\arguments{

\item{theta}{ESF model parameter .} 
\item{n}{Total number of observed elements.}

}

\examples{
x <- ewens.probability_new_next(theta=20, n=100000)
}

