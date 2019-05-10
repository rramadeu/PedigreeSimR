# PedigreeSimR
Wrap-up R package to coordinate PedigreeSim simulations 

```R
devtools::install_github("rramadeu/PedigreeSimR", auth_token = "86fd0062bbfb8632ba148de55eac260f1095553d")
library(PedigreeSimR)
map = 1:100
haplotypes = fake_haplo(100,100,seed=1234)
pedigree = diallel_pedigree(parents=7,popsize=1000,selfs=2)
pedigreesimR(map,haplotypes,pedigree)
```
