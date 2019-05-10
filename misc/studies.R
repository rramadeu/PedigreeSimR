# Studies

grid=expand.grid(totalparents=2:8,
                 totalselfs=0:8,
                 popsize=1000)

grid = grid[c(grid$totalselfs<=grid$totalparents),]

for(i in 1:nrow(grid)){
  x=diallel_pedigree(parents=grid$totalparents[i],popsize=1000,subpopsize=NULL,selfs=grid$totalselfs[i],nextinteger = FALSE)
  grid$popsize[i]=nrow(x)-grid$totalparents[i]
}



map = 1:100
haplotypes = fake_haplo(100,100,seed=1234)
pedigree = diallel_pedigree(parents=7,popsize=1000,selfs=2)
