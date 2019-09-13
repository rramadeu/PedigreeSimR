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
#' @param mapfunction "HALDANE" or "KOSAMBI"
#' @param chromosome string
#' @param sampleHap if TRUE sample haplotypes, pick them in sequence
#' @param seed integer to be used to sample haplotypes
#' @param allownochiasmata numeric
#' @param naturalpairing numeric
#' @param parallelquadrivalents numeric
#' @param pairedcentromeres numeric
#' @param mapwidthpad numeric length of marker code
#' @param GBS if TRUE simulate GBS data and do SNP calling with updog
#' @param GBSavgdepth average depth to sample total number of reads from Poisson distribution
#' @param GBSseq the sequencing error rate (rflexdog inner function)
#' @param GBSod the overdispersion parameter (rflexdog inner function).
#' @param GBSbias the bias parameter  Pr(a read after selected) / Pr(A read after selected) (rflexdog inner function). (rflexdog inner function)
#' @param GBSsnpcall if TRUE performs SNP calling using updog
#' @param GBSnc number of cores for the parallelization for SNP calling
#' @param monoFilter if TRUE filter monomorphic markers from haplotypes
#'
#' @return nothing
#'
#' @examples
#' map = 1:100
#' haplotypes = fake_haplo(m=50,n=100,seed=1234)
#' pedigree = diallel_pedigree(parents=7,popsize=1000,selfs=2)
#' pedigreesimR(map,haplotypes,pedigree)
#'
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
                         sampleHap=TRUE,
                         seed=NULL,
                         allownochiasmata=1,
                         naturalpairing=1,
                         parallelquadrivalents=0,
                         pairedcentromeres=0,
                         mapwidthpad=4,
                         GBS = FALSE,
                         GBSavgdepth=60,
                         GBSsnpcall = FALSE,
                         GBSseq = 0.001,
                         GBSbias = 0.7,
                         GBSod = 0.005,
                         GBSnc = 1,
                         monoFilter = TRUE){

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
    haplotypes = haplotypes[-tmp,]
    mapdf = mapdf[-tmp,]
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
                                       "ALLOWCHIASMATA =",
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
  if(GBS){
    cat("\n Sampling GBS data")

    if(!is.null(seed)) (set.seed(seed))
    sizemat = matrix(rpois(prod(dim(truegenos)),GBSavgdepth),nrow(truegenos),ncol(truegenos))

    refmat = truegenos*0
    for(i in 1:nrow(truegenos)){
      refmat[i,] <- rflexdog(sizevec = as.numeric(sizemat[i,]),
                             geno    = as.numeric(truegenos[i,]),
                             ploidy  = ploidy,
                             seq     = GBSseq,
                             bias    = GBSbias,
                             od      = GBSod)
    }

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

  indnames <- colnames(truegenos)
  marknames <- rownames(truegenos)
  cat("\n Writing PolyOrigin Files")
  ## Formating to PolyOrigin Genotypic Format
  truegenos=cbind(mapdf,read.table(paste0(workingfolder,"/",filename,"pedsim_out_alleledose.dat"),header=TRUE)[,-1])
  write.table(truegenos,file=paste0(workingfolder,"/",filename,"polyorigin_geno.csv"),row.names = FALSE,quote = FALSE,sep=" , ")

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
