# PedigreeSimR
Wrap-up R package to coordinate PedigreeSim simulations ``` 

# Basic usage
```R
library(PedigreeSimR)
map = 1:100
haplotypes = fake_haplo(n=50,m=100,seed=1234)

## Diallel pedigree with 3 parents, no selfs, and total pop size = 90
# Creates 3 fullsib pops from 3 parents (A, B, C), AxB (n=30), AxC (n=30), BxC (n=30)
pedigree = diallel_pedigree(parents=3,popsize=100)

## Diallel pedigree with 3 parents, no selfs, and total pop size = 100
# Creates 3 fullsib pops from 3 parents (A, B, C), AxB (n=334), AxC (n=333), BxC (n=333)
# If the number of populations is not a factor of the popsize
# it creates pops with difference of 1 unit to fit the popsize specified
pedigree = diallel_pedigree(parents=3,popsize=100)

# however, if nextinteger=TRUE, it look for the next integer based on the number of pops
# creates 3 fullsib pops from 3 parents (A, B, C), AxB (n=334), AxC (n=334), BxC (n=334)
pedigree = diallel_pedigree(parents=3,popsize=100,nextinteger = TRUE)

# you can also set the total size of each cross:
pedigree = diallel_pedigree(parents=3,subpopsize=333)

## Diallel pedigree with 6 parents, no selfs, and total pop size = 100
pedigree = diallel_pedigree(parents=6,popsize=1000)

## Diallel pedigree with 6 parents, 2 selfs, and total pop size = 100
pedigree = diallel_pedigree(parents=6,selfs=2,popsize=100)

## Diallel pedigree with 6 parents and 2 selfs, and each pop (cross) with size = 100
pedigree = diallel_pedigree(parents=6,selfs=2,subpopsize=100)

## Single-Round-Robin with 6 parents and total pop size = 100
pedigree = round_pedigree(parents=6,popsize=100)

## Single-Round-Robin with 6 parents and total pop size = 100
pedigree = round_pedigree(parents=6,popsize=100,nextinteger=TRUE)

## Single-Round-Robin with 5 parents and each pop (cross) with size = 100
pedigree = round_pedigree(parents=6,subpopsize=100)

## Simulating genotypes
pedigreesimR(map,haplotypes,pedigree,ploidy=4)

## Simulating genotypes and GBS data (avg depth=60, seqerror = 0.001, allelic bias = 0.7, overdispersion = 0.005)
pedigreesimR(map,haplotypes,pedigree,ploidy=4,GBS=TRUE,GBSavgdepth = 60,GBSseq = 0.001,GBSbias = 0.7,GBSod = 0.005)

## Simulating genotypes with GBS data and doing SNP calling with updog with 2 cores considering the general model (it takes a while)
pedigreesimR(map,haplotypes,pedigree,ploidy=4,GBS=TRUE,GBSavgdepth = 60,GBSseq = 0.001,GBSbias = 0.7,GBSod = 0.005,GBSsnpcall=TRUE,GBSnc = 2)
```
