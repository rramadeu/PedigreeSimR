#' PedigreeSim Management
#'
#' Wrap-up function to run PedigreeSim software, to simulate GBS data, and to do SNP calling with updog package
#'
#' @param map vector of the marker position.
#' @param haplotypes matrix with haplotype information to be draw
#' @param pedigree pedigree data frame
#' @param centromere numeric with its position
#' @param prefPairing numeric
#' @param quadrivalents numeric
#' @param ploidy integer even number
#' @param workingfolder string with folder name to write input/output files
#' @param filename string with filename id to put in the beggining of the created files
#' @param mapfunction "HALDANE" or "KOSAMBI"
#' @param chromosome string
#' @param sampleHap if TRUE sample haplotypes, if FALSE pick them in sequence
#' @param seed integer to be used to sample haplotypes
#' @param allownochiasmata numeric
#' @param naturalpairing numeric
#' @param parallelquadrivalents numeric
#' @param pairedcentromeres numeric
#' @param mapwidthpad numeric length of marker code
#' @param epsilon vector with two epsilon values for error for genotypic assigment, vector with two values, the first is the parents' epsilon, the second is the offsprings
#' @param missingFreq vector with two missingFreq values to sample missing values for the genotypes, the first is the parents' frequency, the second is the offsprings
#' @param GBS if TRUE simulate GBS data and do SNP calling with updog
#' @param GBSavgdepth average depth to sample total number of reads from Poisson distribution
#' @param GBSseq the sequencing error rate (rflexdog inner function)
#' @param GBSod the overdispersion parameter (rflexdog inner function).
#' @param GBSbias the bias parameter  Pr(a read after selected) / Pr(A read after selected) (rflexdog inner function). (rflexdog inner function)
#' @param GBSsnpcall if TRUE performs SNP calling using updog
#' @param GBSnc number of cores for the parallelization for SNP calling
#' @param monoFilter if TRUE filter monomorphic markers from haplotypes
#' @param trackErrorSim if TRUE create a spreadsheet with the simulate error positions given epsilon and missingData (1 if error, 0 if not)
#' @param justEssential it will just run the simulation and avoid GBS/Dosage Calling/Parental haplotype
#'
#' @return nothing
#'
#' @examples
#' \dontrun{
#' map = 1:100
#' haplotypes = fake_haplo(m=50,n=100,seed=1234)
#' pedigree = diallel_pedigree(parents=7,popsize=1000,selfs=2)
#' pedigreesimR(map,haplotypes,pedigree)
#' }
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export
#'


pedigreesimR <- function(map,
                         haplotypes,
                         pedigree,
                         centromere=NULL,
                         prefPairing=0,
                         quadrivalents=0,
                         ploidy=4,
                         workingfolder="PedigreeSimR_files",
                         filename = "",
                         mapfunction="HALDANE",
                         chromosome="A",
                         sampleHap=FALSE,
                         seed=NULL,
                         allownochiasmata=1,
                         naturalpairing=0,
                         parallelquadrivalents=0,
                         pairedcentromeres=0,
                         mapwidthpad=4,
                         epsilon=c(0,0),
                         missingFreq=c(0,0),
                         GBS = FALSE,
                         GBSavgdepth=60,
                         GBSsnpcall = FALSE,
                         GBSseq = 0.001,
                         GBSbias = 0.7,
                         GBSod = 0.005,
                         GBSnc = 1,
                         monoFilter = TRUE,
                         trackErrorSim = FALSE,
                         justEssential = FALSE){

  ## Creating map file
  mapdf = data.frame(marker=paste0(chromosome,"_",str_pad(1:length(map),width = mapwidthpad,side = "left",pad = "0")),
                     chromosome=chromosome,
                     position=map)

  ## Creating chromosome file
  if(is.null(centromere))
    centromere=max(map)/2

  chrdf = data.frame(chromosome=chromosome,
                     length=max(map),
                     centromere=centromere,
                     prefPairing=prefPairing,
                     quadrivalents=quadrivalents)

  ## Creating founder file
  founders = pedigree[pedigree[,2]=="NA",1]
  totalfounders = length(founders)
  if(!is.null(seed)) (set.seed(seed))
    sampledhaplos = sample(1:ncol(haplotypes),totalfounders*ploidy)
  if(sampleHap){
    haplotypes = haplotypes[,sampledhaplos]
    }else{
    haplotypes = haplotypes[,1:c(totalfounders*ploidy)]
  }

  if(monoFilter){
    tmp = which(apply(haplotypes,1,var) == 0)
    if(length(tmp)>0){
      haplotypes = haplotypes[-tmp,]
      mapdf = mapdf[-tmp,]
    }
  }

  colnames(haplotypes) = paste0(rep(founders,each=ploidy),"_",1:ploidy)
  founderdf = data.frame(marker=mapdf$marker,
                         haplotypes)

  ## Parameter file
  parameterdf = data.frame(parameter=c("PLOIDY =",
                                       "MAPFUNCTION =",
                                       "MISSING =",
                                       "CHROMFILE =",
                                       "PEDFILE =",
                                       "MAPFILE =",
                                       "FOUNDERFILE =",
                                       "OUTPUT =",
                                       "ALLOWNOCHIASMATA =",
                                       "NATURALPAIRING =",
                                       "PARALLELQUADRIVALENTS =",
                                       "PAIREDCENTROMERES ="),
                           values=c(ploidy,
                                    mapfunction,
                                    "NA",
                                    paste0(workingfolder,"/",filename,"pedsim_input.chrom"),
                                    paste0(workingfolder,"/",filename,"pedsim_input.ped"),
                                    paste0(workingfolder,"/",filename,"pedsim_input.map"),
                                    paste0(workingfolder,"/",filename,"pedsim_input.founder"),
                                    paste0(workingfolder,"/",filename,"pedsim_out"),
                                    allownochiasmata,
                                    naturalpairing,
                                    parallelquadrivalents,
                                    pairedcentromeres)
  )

  ## Writing files
  if(workingfolder != getwd())
    if(is.na(match(workingfolder,list.files())))
      dir.create(workingfolder)

  write.table(chrdf,file=paste0(workingfolder,"/",filename,"pedsim_input.chrom"),row.names = FALSE,quote = FALSE)
  write.table(pedigree,file=paste0(workingfolder,"/",filename,"pedsim_input.ped"),row.names = FALSE,quote = FALSE)
  write.table(mapdf,file=paste0(workingfolder,"/",filename,"pedsim_input.map"),row.names = FALSE,quote = FALSE)
  write.table(founderdf,file=paste0(workingfolder,"/",filename,"pedsim_input.founder"),row.names = FALSE,quote = FALSE)
  write.table(parameterdf,file=paste0(workingfolder,"/",filename,"pedsim_input.par"),row.names = FALSE, col.names=FALSE,quote = FALSE)

  cat("\n Simulating with PedigreeSim \n")
  system(
    paste0("java -jar \"",system.file(package = "PedigreeSimR"),"/PedigreeSim2/PedigreeSim.jar\" ",workingfolder,"/",filename,"pedsim_input.par")
  )

  ## Formating to PolyOrigin Pedigree
  levels(pedigree$Parent1)[which(levels(pedigree$Parent1)=="NA")] = "0"
  levels(pedigree$Parent2)[which(levels(pedigree$Parent2)=="NA")] = "0"
  pedigree$Population = substr(pedigree$Name,1,3)
  pedigree$Population[which(nchar(pedigree$Population)==1)] = 0
  pedigree$Population = as.numeric(as.factor(pedigree$Population))-1
  pedigree$Ploidy=ploidy
  pedigree=pedigree[,c(1,4,2,3,5)]
  names(pedigree) = c("Individual","Population","MotherID","FatherID","Ploidy")
  write.table(pedigree,file=paste0(workingfolder,"/",filename,"polyorigin_pedigree.csv"),row.names = FALSE,quote = FALSE,sep=" , ")
  truegenos = read.table(paste0(workingfolder,"/",filename,"pedsim_out_alleledose.dat"),header=TRUE)[,-1]


  ## Sampling sequencing depth
  ## Adding genotyping errors
  total.parents = length(founders)

  if(sum(epsilon>0)){
      truegenos.par = as.matrix(truegenos[,1:total.parents])
      truegenos.off = as.matrix(truegenos[,-c(1:total.parents)])
      colnames(truegenos.par) = founders

      E.par = matrix(rbinom(prod(c(length(truegenos.par),1)),1,epsilon[1]),nrow=length(truegenos.par))
      E.par = matrix(rbinom(prod(dim(truegenos.par)),1,epsilon[1]),nrow=nrow(truegenos.par))

      if(sum(E.par)>0){
        which.E.par = which(E.par==1,arr.ind = TRUE)
        for(i in 1:nrow(which.E.par))
          truegenos.par[which.E.par[i,1],which.E.par[i,2]] <- sample(c(0:ploidy)[-(truegenos.par[which.E.par[i,1],which.E.par[i,2]]+1)],1)
      }
      E.off = matrix(rbinom(prod(dim(truegenos.off)),1,epsilon[2]),nrow=nrow(truegenos.off))
      if(sum(E.off)>0){
        which.E.off = which(E.off==1,arr.ind = TRUE)
        for(i in 1:nrow(which.E.off))
          truegenos.off[which.E.off[i,1],which.E.off[i,2]] <- sample(c(0:ploidy)[-(truegenos.off[which.E.off[i,1],which.E.off[i,2]]+1)],1)
      }
      truegenos.eps = cbind(truegenos.par,truegenos.off)
  }


  if(sum(missingFreq>0)){
    truegenos.par = as.matrix(truegenos[,1:total.parents])
    truegenos.off = as.matrix(truegenos[,-c(1:total.parents)])
    colnames(truegenos.par) = founders

    F.par = matrix(rbinom(prod(dim(truegenos.par)),1,missingFreq[1]),nrow=nrow(truegenos.par))
    if(sum(F.par)>0)
      truegenos.par[which(F.par==1,arr.ind=1)] = NA

    F.off = matrix(rbinom(prod(dim(truegenos.off)),1,missingFreq[2]),nrow=nrow(truegenos.off))
    if(sum(F.off)>0)
      truegenos.off[which(F.off==1,arr.ind=1)] = NA
    truegenos.NA = cbind(truegenos.par,truegenos.off)
    truegenos.NA = is.na(truegenos.NA)
  }

  if(sum(missingFreq>0)){
    if(sum(epsilon==0))
      truegenos.eps = truegenos
    truegenosNA = truegenos.eps
    truegenosNA[truegenos.NA] = NA
    ## Formating to PolyOrigin Genotypic Format
    write.table(cbind(mapdf,truegenosNA),file=paste0(workingfolder,"/",filename,"polyorigin_geno_snparray.csv"),row.names = FALSE,quote = FALSE,sep=" , ")
    #if(trackErrorSim){
    #  write.table(cbind(mapdf,ifelse(is.na(truegenosNA),1,0)),file=paste0(workingfolder,"/",filename,"polyorigin_geno_missingdata_track.csv"),row.names = FALSE,quote = FALSE,sep=" , ")
    #}
  }else{
    if(sum(epsilon==0))
      truegenos.eps = truegenos
    write.table(cbind(mapdf,truegenos.eps),file=paste0(workingfolder,"/",filename,"polyorigin_geno_snparray.csv"),row.names = FALSE,quote = FALSE,sep=" , ")
  }

  if(GBS){
    cat("\n Sampling GBS data")

    if(!is.null(seed)) (set.seed(seed))
    sizemat = matrix(rpois(prod(dim(truegenos)),GBSavgdepth),nrow(truegenos),ncol(truegenos))

    refmat = truegenos.eps*0
    for(i in 1:nrow(truegenos)){
      refmat[i,] <- rflexdog(sizevec = as.numeric(sizemat[i,]),
                             geno    = as.numeric(truegenos[i,]),
                             ploidy  = ploidy,
                             seq     = GBSseq,
                             bias    = GBSbias,
                             od      = GBSod)
    }

    if(sum(missingFreq>0))
      refmat[truegenos.NA] = NA

    if(GBSsnpcall){
      if(GBSnc==1){
        cat("\n Doing SNP calling")
        geno = NULL
        for(i in 1:nrow(truegenos)){
          fout <- flexdog(as.numeric(refmat[i,]),
                             sizemat[i,],
                             ploidy=ploidy,
                             verbose=FALSE,
                     model="norm")
          geno <- cbind(geno,cbind(fout$geno,apply(round(fout$postmat,3),1,paste0,collapse="|")))
        }
      }else{
        cat(paste("\n Doing SNP calling with",GBSnc,"cores"))
        ## number of cores.
        ## You should change this for your specific computing environment.
        cl <- parallel::makeCluster(GBSnc)
        ngenes <- nrow(refmat)
        doParallel::registerDoParallel(cl = cl)
        stopifnot(foreach::getDoParWorkers() > 1) ## make sure cluster is set up.
        geno <- foreach(i = 1:ngenes,
                        .combine = cbind,
                        .export = "flexdog") %dopar% {
                          ## fit flexdog
                          fout <- flexdog(refvec  = as.numeric(refmat[i,]),
                                          sizevec = sizemat[i,],
                                          ploidy  = ploidy,
                                          model   = "norm",
                                          verbose = FALSE)
                          cbind(fout$geno,apply(round(fout$postmat,3),1,paste0,collapse="|"))
                        }

        stopCluster(cl)
      }
    }
  }

  cat("\n Writing PolyOrigin Files")

  indnames <- colnames(truegenos)
  marknames <- rownames(truegenos)

  ## Formating to PolyOrigin Genotypic Format
  truegenos=cbind(mapdf,read.table(paste0(workingfolder,"/",filename,"pedsim_out_alleledose.dat"),header=TRUE)[,-1])
  write.table(truegenos,file=paste0(workingfolder,"/",filename,"polyorigin_geno.csv"),row.names = FALSE,quote = FALSE,sep=" , ")

  if(justEssential){
    cat("\n Done!")
    return()
  }


  ## Extracting true haplotype information given pedigree
  truehaplos = read.table(paste0(workingfolder,"/",filename,"pedsim_out_founderalleles.dat"),header=TRUE)[,-1]
  #total.parents = ncol(haplotypes)/ploidy
  truehaplos = truehaplos[,-(1:(ploidy*total.parents))]

  ## Haplo within info
  truehaploW = (truehaplos %% ploidy) + 1

  ## Trick to track when self or not and multiple by 4 the correct parent
  track.selfs = as.character(pedigree$MotherID)!=as.character(pedigree$FatherID)
  track.selfs = track.selfs[-(1:total.parents)]
  track.selfs = c(sapply(seq_along(track.selfs), function(i) append(rep(FALSE,length(track.selfs))[i], track.selfs[i], i)))
  track.selfs = track.selfs*ploidy
  track.selfs = rep(track.selfs,each=ploidy/2)
  track.selfs = matrix(rep(track.selfs,nrow(truehaploW)),nrow=nrow(truehaploW),byrow = TRUE)

  ## truehaploW back
  truehaploW = truehaploW + track.selfs
  inds <- as.character(pedigree$Individual)[-c(1:total.parents)]
  inds.hap <- rep(inds,each=ploidy)

  truehaploCollapsed = NULL
  for(i in 1:length(inds)){
    ind.match = which((match(inds.hap,inds[i]))==1)
    truehaploCollapsed = cbind(truehaploCollapsed,apply(truehaploW[,ind.match],1,paste,collapse="|"))
  }
  colnames(truehaploCollapsed) = inds

  ## Extracting geno position from parents haplotype
  truehaplos = read.table(paste0(workingfolder,"/",filename,"pedsim_input.founder"),header=TRUE)[,-1]
  truehaplos = truehaplos+1
  parents = as.character(pedigree$Individual[1:total.parents])
  parents.hap <- rep(parents,each=ploidy)
  parenthaploCollapsed = NULL
  for(i in 1:length(parents)){
    ind.match = which((match(parents.hap,parents[i]))==1)
    parenthaploCollapsed = cbind(parenthaploCollapsed,apply(truehaplos[,ind.match],1,paste,collapse="|"))
  }
  colnames(parenthaploCollapsed) = parents
  haploexport = cbind(truegenos[,1:3],parenthaploCollapsed,truehaploCollapsed)
  write.table(haploexport,file=paste0(workingfolder,"/",filename,"polyorigin_truevalue_ancestral.csv"),row.names = FALSE,quote = FALSE,sep=" , ")

  ## Extracting geno position from parents haplotype (just presence/absence information)
  inds <- as.character(pedigree$Individual)
  inds.hap <- rep(inds,each=ploidy)

  truehaplos = read.table(paste0(workingfolder,"/",filename,"pedsim_out_founderalleles.dat"),header=TRUE)[,-1]
  genofounders = read.table(paste0(workingfolder,"/",filename,"pedsim_input.founder"),header=TRUE)[,-1]
  genofounders = genofounders+1
  truehaplos10 = truehaplos
  for(i in 1:nrow(truehaplos)){
    tmp = factor(truehaplos[i,],levels=as.character(0:(ploidy*total.parents-1)))
    levels(tmp) = as.numeric(genofounders[i,])
    truehaplos10[i,] = tmp
  }

  truehaploCollapsed = NULL
  for(i in 1:length(inds)){
    ind.match = which((match(inds.hap,inds[i]))==1)
    truehaploCollapsed = cbind(truehaploCollapsed,apply(truehaplos10[,ind.match],1,paste,collapse="|"))
  }
  colnames(truehaploCollapsed) = inds
  haploexport = cbind(truegenos[,1:3],truehaploCollapsed)
  write.table(haploexport,file=paste0(workingfolder,"/",filename,"polyorigin_truevalue.csv"),row.names = FALSE,quote = FALSE,sep=" , ")

  if(sum(epsilon>0)){
    if(trackErrorSim){
      haploexport = haploexport[,-c(1:3)]
      truegenos.track = truegenos[,-c(1:3)]
      truegenos.eps.track=ifelse(truegenos.eps!=truegenos.track,1,0)
      track.index=which(truegenos.eps.track==1)
      track.index.arr=which(truegenos.eps.track==1,arr.ind = TRUE)
      truegenos.eps.track=data.frame(mapdf[track.index.arr[,1],],
                                     pedigree[track.index.arr[,2],c(1,3,4)],
                                     wrongdosage=as.matrix(truegenos.eps)[track.index],
                                     truedosage=as.matrix(truegenos.track)[track.index],
                                     truephases=as.matrix(haploexport)[track.index],
                                     na=as.matrix(truegenos.NA)[track.index])
      truegenos.eps.track = truegenos.eps.track[-truegenos.eps.track$na,-10]
      write.table(truegenos.eps.track,file=paste0(workingfolder,"/",filename,"polyorigin_geno_snparray_errortrack.csv"),row.names = FALSE,quote = FALSE,sep=" , ")
    }
  }

  ## With SNPCalling/Geno Error
  if(GBS){

    altmat=as.matrix(sizemat-refmat)
    refmat=as.matrix(refmat)
    count=paste0(altmat,"|",refmat)
    dim(count) <- dim(altmat)
    colnames(count) <- colnames(altmat)
    rownames(count) <- rownames(altmat)
    callinggenos=cbind(truegenos[,1:3],count)
    write.table(callinggenos,file=paste0(workingfolder,"/",filename,"polyorigin_geno_GBS.csv"),row.names = FALSE,quote = FALSE,sep=" , ")
  }
    if(GBSsnpcall){
        postmat <- t(geno[,seq(from=2,to=ncol(geno),by=2)])
        geno <- t(geno[,seq(from=1,to=ncol(geno),by=2)])
        colnames(postmat) <- colnames(geno) <- indnames
        callinggenos=cbind(truegenos[,1:3],geno)
        callingpostmat=cbind(truegenos[,1:3],postmat)
        write.table(callinggenos,file=paste0(workingfolder,"/",filename,"polyorigin_geno_GBS_updog.csv"),row.names = FALSE,quote = FALSE,sep=" , ")
       write.table(callingpostmat,file=paste0(workingfolder,"/",filename,"polyorigin_geno_GBS_updog_postmat.csv"),row.names = FALSE,quote = FALSE,sep=" , ")
  }
    cat("\n Done!")
}
