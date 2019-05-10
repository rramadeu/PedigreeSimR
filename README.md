# PedigreeSimR
Wrap-up R package to coordinate PedigreeSim simulations 

# Installing
```R
devtools::install_github("rramadeu/PedigreeSimR", 
                         auth_token = "86fd0062bbfb8632ba148de55eac260f1095553d")
``` 

# Basic usage
```R
library(PedigreeSimR)
map = 1:100
haplotypes = fake_haplo(n=50,m=100,seed=1234)

## Diallel pedigree with 7 parents and 2 selfs
diallel7 = diallel_pedigree(parents=7,popsize=1000,selfs=2)
pedigreesimR(map,haplotypes,diallel7,ploidy=4)

## Single Round Robin pedigree with 7 parents
round7 = round_pedigree(parents=7,popsize=1000)
pedigreesimR(map,haplotypes,round7,ploidy=4)
```
