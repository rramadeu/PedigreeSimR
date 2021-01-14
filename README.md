# PedigreeSimR
Wrapper R package to coordinate `PedigreeSim` simulations for classic breeder mating designs in outcrossing plants. It also simulates genotypic data from GBS. This package uses some `updog` package functions for the simulation of the genotypes and GBS data.

# Citation

I created PedigreeSimR as an auxiliary package to do the simulations for these two manuscripts (pre-prints):

Amadeu, R. R., Munoz, P., Zheng, C., & Endelman, J. B. (2020). QTL Mapping in Outbred Tetraploid (and Diploid) Diallel Populations. bioRxiv. [https://doi.org/10.1101/2020.12.18.423479](https://doi.org/10.1101/2020.12.18.423479)

Zheng, C., Amadeu, R. R., Munoz, P., & Endelman, J. B. (2020). Haplotype reconstruction in connected tetraploid F1 populations. bioRxiv. [https://doi.org/10.1101/2020.12.18.423519](https://doi.org/10.1101/2020.12.18.423519)

As there is no proper manuscript, you can cite either of those two pre-prints as reference source by now.

# References:

[PedigreeSim](https://www.wur.nl/en/show/Software-PedigreeSim.htm):

Voorrips, R. E., & Maliepaard, C. A. (2012). The simulation of meiosis in diploid and tetraploid organisms using various genetic models. BMC bioinformatics, 13(1), 248.

[updog](https://CRAN.R-project.org/package=updog):

Gerard D., Ferrão L., Garcia A., Stephens M. (2018). Genotyping Polyploids from Messy Sequencing Data. Genetics, 210(3), 789–807. ISSN 0016-6731

Gerard D., & Ferrão L. (2020). Priors for Genotyping Polyploids. Bioinformatics, 36(6), 1795–1800. ISSN 1367-4803



# Basic usage in R
```R
## Installing this PedigreeSimR package
devtools::install_github("rramadeu/PedigreeSimR")
library(PedigreeSimR)

## Creating a fake map (here you add a real map with haplotypes)
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
