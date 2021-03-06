% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pedigreesimR.R
\name{pedigreesimR}
\alias{pedigreesimR}
\title{PedigreeSim Management}
\usage{
pedigreesimR(
  map,
  haplotypes,
  pedigree,
  centromere = NULL,
  prefPairing = 0,
  quadrivalents = 0,
  ploidy = 4,
  workingfolder = "PedigreeSimR_files",
  filename = "",
  mapfunction = "HALDANE",
  chromosome = "A",
  sampleHap = FALSE,
  seed = NULL,
  allownochiasmata = 1,
  naturalpairing = 0,
  parallelquadrivalents = 0,
  pairedcentromeres = 0,
  mapwidthpad = 4,
  epsilon = c(0, 0),
  missingFreq = c(0, 0),
  GBS = FALSE,
  GBSavgdepth = 60,
  GBSsnpcall = FALSE,
  GBSseq = 0.001,
  GBSbias = 0.7,
  GBSod = 0.005,
  GBSnc = 1,
  monoFilter = TRUE,
  trackErrorSim = FALSE,
  justEssential = FALSE
)
}
\arguments{
\item{map}{vector of the marker position.}

\item{haplotypes}{matrix with haplotype information to be draw}

\item{pedigree}{pedigree data frame}

\item{centromere}{numeric with its position}

\item{prefPairing}{numeric}

\item{quadrivalents}{numeric}

\item{ploidy}{integer even number}

\item{workingfolder}{string with folder name to write input/output files}

\item{filename}{string with filename id to put in the beggining of the created files}

\item{mapfunction}{"HALDANE" or "KOSAMBI"}

\item{chromosome}{string}

\item{sampleHap}{if TRUE sample haplotypes, if FALSE pick them in sequence}

\item{seed}{integer to be used to sample haplotypes}

\item{allownochiasmata}{numeric}

\item{naturalpairing}{numeric}

\item{parallelquadrivalents}{numeric}

\item{pairedcentromeres}{numeric}

\item{mapwidthpad}{numeric length of marker code}

\item{epsilon}{vector with two epsilon values for error for genotypic assigment, vector with two values, the first is the parents' epsilon, the second is the offsprings}

\item{missingFreq}{vector with two missingFreq values to sample missing values for the genotypes, the first is the parents' frequency, the second is the offsprings}

\item{GBS}{if TRUE simulate GBS data and do SNP calling with updog}

\item{GBSavgdepth}{average depth to sample total number of reads from Poisson distribution}

\item{GBSsnpcall}{if TRUE performs SNP calling using updog}

\item{GBSseq}{the sequencing error rate (rflexdog inner function)}

\item{GBSbias}{the bias parameter  Pr(a read after selected) / Pr(A read after selected) (rflexdog inner function). (rflexdog inner function)}

\item{GBSod}{the overdispersion parameter (rflexdog inner function).}

\item{GBSnc}{number of cores for the parallelization for SNP calling}

\item{monoFilter}{if TRUE filter monomorphic markers from haplotypes}

\item{trackErrorSim}{if TRUE create a spreadsheet with the simulate error positions given epsilon and missingData (1 if error, 0 if not)}

\item{justEssential}{it will just run the simulation and avoid GBS/Dosage Calling/Parental haplotype}
}
\value{
nothing
}
\description{
Wrap-up function to run PedigreeSim software, to simulate GBS data, and to do SNP calling with updog package
}
\examples{
\dontrun{
map = 1:100
haplotypes = fake_haplo(m=50,n=100,seed=1234)
pedigree = diallel_pedigree(parents=7,popsize=1000,selfs=2)
pedigreesimR(map,haplotypes,pedigree)
}
}
\author{
Rodrigo R Amadeu, \email{rramadeu@gmail.com}
}
