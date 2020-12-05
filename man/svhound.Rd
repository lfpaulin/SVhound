\name{svhound}
\docType{package}
\alias{svhound}
\alias{svhound-package}
\title{SVhound}

\description{ TODO }


\usage{ sv_hound(structuralVariantsData="data/my_sv_dataset.svalleles", window_size=10000, output_prefix="my_sv_dataset") }

\arguments{

    \item{structuralVariantsData}{Input dataset, it can be a VCF file for which parsing instructions will be given or a TSV file output from vcf_parser_for_svhound.py}
    \item{window_size}{Integer that denotes the size of the windows used to discretize the genome.}
    \item{output_prefix}{Prefix of the output file.}
    \item{subsample}{Number of genomes/individuals to be sampled from the 'structuralVariantsData'. It needs to be an integer greater than zero. Default=NULL}
    \item{usePSF}{Parameter to use the Pitman Sampling Formula instead of the Ewens Sampling Formula as the model. TURE/FALSE, default=FALSE}
    \item{giveExampleData}{Returns an example of the TSV file output from vcf_parser_for_svhound.py loaded in R. Default=FALSE}
    \item{runExample}{Runs the example dataset. Default=FALSE}

}

\keyword{package}

