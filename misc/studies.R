# Studies

setwd("~/Documents/PedigreeSimR_misc/")
map = seq(1,3,0.5)
haplotypes = fake_haplo(m=length(map),n=100,seed=1234)
pedigree = diallel_pedigree(parents=7,popsize=1000,selfs=2)
pedigreesimR(map,haplotypes,pedigree,ploidy=4,GBS=FALSE,GBSnc = 2)
pedigreesimR(map,haplotypes,pedigree,ploidy=4,GBS=TRUE,GBSnc = 2)
pedigreesimR(map,haplotypes,pedigree,ploidy=4,GBS=TRUE,GBSsnpcall=TRUE,GBSnc = 1)
pedigreesimR(map,haplotypes,pedigree,ploidy=4,GBS=TRUE,GBSsnpcall=TRUE,GBSnc = 2)

## TetraOrigin Example
setwd("~/Documents/PedigreeSimR_misc/")
map = seq(1,3,0.5)
haplotypes = fake_haplo(m=length(map),n=100,seed=1234)
pedigree = diallel_pedigree(parents=7,popsize=1000,selfs=2)
pedigreesimR(map,haplotypes,pedigree,ploidy=4,GBS=TRUE,GBSsnpcall = FALSE, GBSnc = 2)


library(updog)
map=map
haplotypes=haplotypes
pedigree=pedigree
centromere=NULL
prefPairing=0
quadrivalents=0
ploidy=2
workingfolder="PedigreeSimR_files"
mapfunction="HALDANE"
chromosome="A"
seed=NULL
allownochiasmata=1
naturalpairing=1
parallelquadrivalents=0
pairedcentromeres=0
mapwidthpad=4
GBS = TRUE
GBSseq = 0.001
GBSbias = 0.7
GBSod = 0.005
GBSnc = 2
GBSsnpcall=FALSE

set.seed(1)
library(updog)
nind    <- 100
ploidy  <- 6
sizevec <- round(stats::runif(n   = nind,
                              min = 50,
                              max = 200))
true_geno <- rgeno(n      = nind,
                   ploidy = ploidy,
                   model  = "f1",
                   p1geno = 4,
                   p2geno = 5)

refvec <- rflexdog(sizevec = sizevec,
                   geno    = true_geno,
                   ploidy  = ploidy,
                   seq     = 0.001,
                   bias    = 0.7,
                   od      = 0.005)

plot_geno(refvec  = refvec,
          sizevec = sizevec,
          ploidy  = ploidy,
          bias    = 0.7,
          seq     = 0.001,
          geno    = true_geno)


