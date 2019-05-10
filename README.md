# PedigreeSimR
Wrap-up R package to coordinate PedigreeSim simulations 


map = 1:100
haplotypes = fake_haplo(100,100,seed=1234)
pedigree = diallel_pedigree(parents=7,popsize=1000,selfs=2)
pedigreesimR(map,haplotypes,pedigree)
