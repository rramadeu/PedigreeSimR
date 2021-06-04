#' Simulate QTL + residual effects
#'
#' Simulates phenotypic data on top of the pedigreeSimR output and save then as diaQTL csv format
#'
#' @param parents number of parents in the pedigree
#' @param ploidy integer even number
#' @param workingfolder string with folder name to load the input and write output files
#' @param QTLmarker name of marker to be used as QTL, if NULL a marker is sampled randomly (default)
#' @param QTLh2 numeric between 0 and 1 with the QTL heritability
#' @param run_diaQTL logical, if TRUE runs diaQTL scan1 and fitQTL functions in the data.
#'
#' @return if run_diaQTL=TRUE, return a list with diaQTL results. Otherwise, returns NULL.
#'
#' @examples
#' \dontrun{
#' map = 1:100
#' haplotypes = fake_haplo(m=50,n=100,seed=1234)
#' pedigree = diallel_pedigree(parents=3,popsize=100)
#' pedigreesimR(map,haplotypes,pedigree,workingfolder="PedigreeSimR_files")
#' QTLsim(parents=3,
#'        ploidy=4,
#'        workingfolder="PedigreeSimR_files",
#'        QTLmarker=NULL,
#'        QTLh2=0,
#'        run_diaQTL=FALSE)
#' }
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#' @export
#'

QTLsim <- function(parents=3,
                   ploidy=4,
                   workingfolder="PedigreeSimR_files",
                   QTLmarker=NULL,
                   QTLh2=0.3,
                   run_diaQTL=FALSE){
  init.wd = getwd()
  time_init = Sys.time()
  ## Getting the states notation
  #F1states = diaQTL::F1codes$State
  F1states1 = stringr::str_split(F1codes[,2],"-",simplify=TRUE)
  F1states2 = F1states1[,c(2,1,3,4)]
  F1states3 = F1states1[,c(1,2,4,3)]
  F1states4 = F1states1[,c(2,1,4,3)]
  F1states = rbind(F1states1,F1states2,F1states3,F1states4)
  F1states = apply(F1states,1,paste0,collapse="|")
  F1code = rep(1:100,4)

  pedigree=read.table(paste0(workingfolder,"/pedsim_input.ped"),header = TRUE)
  colnames(pedigree) <- c("id","mother","father")
  founders <- read.csv(paste0(workingfolder,"/polyorigin_truevalue_ancestral.csv"),colClasses="character",sep=",")
  rownames(founders) <- founders[,1]
  foundersmap <- founders[,1:3]
  founders <- founders[,-c(1:(3+parents))]
  founders <- apply(founders,2,function(x) gsub('\\s+', '',x)) #rm white space
  founders <- apply(founders,2,function(x) return(F1code[match(x,F1states)]))
  for(j in 1:ncol(founders))
    founders[,j] <- paste0(founders[,j],"=>1")
  founders <- cbind(foundersmap,founders)
  colnames(founders) <- gsub('\\s+', '',colnames(founders)) #rm white space
  founders[,1] <- gsub('\\s+', '',founders[,1]) #rm white space
  founders[,2] <- gsub('\\s+', '',founders[,2]) #rm white space
  founders[,3] <- gsub('\\s+', '',founders[,3]) #rm white space
  names(founders)[1:3] = c("marker","chrom","cM")
  write.csv(founders,file=paste0(workingfolder,"/geno.csv"),row.names=FALSE, quote=FALSE)

  colnames(pedigree) <- c("id","mother","father")
  write.csv(pedigree[-c(1:parents),],file=paste0(workingfolder,"/QTLped.csv"),row.names=FALSE,quote=FALSE)
  PopSize=nrow(pedigree[-c(1:parents),])
  if(is.null(QTLmarker)){
    sampleQTL = sample(2:(nrow(founders)-1),1)
  }else{
    sampleQTL = match(QTLmarker,founders$markers)
  }

  cat("Genotypic data saved as geno.csv (it contains the QTL marker) \n")
  cat("Pedigree data saved as QTLped.csv \n")

  setwd(paste0(workingfolder))
  data1 <- read_data(genofile="geno.csv",
                     pedfile="QTLped.csv",
                     ploidy=4,
                     dominance=1)
  cat("Simulating the QTL...")
  #do not sample ending markers
  sampleQTLname = data1@map$marker[sampleQTL]
  ## qtlAlleles randomly assigned
  qtlAlleles = rnorm(ploidy*parents,0,10)
  #qtlAlleles = sample(qtlAlleles,ploidy*parents,replace=TRUE)
  qtlEffects = as.numeric(data1@geno[[sampleQTL]][[1]] %*% qtlAlleles)
  if(QTLh2>0){
    for(i in 1:1000){
      residualDummy = rnorm(PopSize,0,10)
      scalingh2= sqrt( (1-QTLh2)/QTLh2 * var(as.numeric(qtlEffects)))

      residualEffects = residualDummy / (sd(residualDummy)/scalingh2)

      pheno = qtlEffects + residualEffects
      pheno = as.vector(pheno)
      h2simul = var(qtlEffects)/var(pheno)
      if(abs(h2simul - QTLh2)<0.0001)
        break
    }
    qtlAllelesAll = qtlAlleles
    qtlEffectsAll = qtlEffects
  }else{
    qtlAllelesAll=qtlEffectsAll=0
    pheno = rnorm(PopSize,0,10)
  }

  write.csv(data.frame(id=rownames(data1@X.GCA),pheno1=pheno),file="QTLpheno.csv",row.names=FALSE,quote=FALSE)
  write.csv(founders[-sampleQTL,],file=paste0("QTLgeno.csv"),row.names=FALSE)
  rm(data1);gc()
  cat("Genotypic data saved as QTLgeno.csv \n")
  cat("Phenotypic data saved as QTLpheno.csv \n")
  if(run_diaQTL){
    data1 <- read_data(genofile="QTLgeno.csv",
                       pedfile="QTLped.csv",
                       phenofile="QTLpheno.csv",
                       ploidy=4,
                       dominance=2)
    par1 <- set_params(data = data1,
                       trait = "pheno1",
                       q=0.5,
                       r=0.1)
    par2 <- set_params(data = data1,
                       trait = "pheno1",
                       q=0.05,
                       r=0.025)
    flankingMarkers = c(data1@map[sampleQTL-1,1],data1@map[sampleQTL,1])
    ans <- scan1(data=data1,trait="pheno1",params=list(burnIn = max(par1[,1]), nIter = max(par1[,2])),n.core=1,dominance = 1)
    qtlfitDown <- fitQTL(data=data1,trait="pheno1",params = list(burnIn = max(par2[,1]), nIter = max(par2[,2])),polygenic = FALSE,
                         qtl=data.frame(marker=flankingMarkers[1],dominance=1))
    qtlfitUp <- fitQTL(data=data1,trait="pheno1",params = list(burnIn = max(par2[,1]), nIter = max(par2[,2])),polygenic = FALSE,
                       qtl=data.frame(marker=flankingMarkers[2],dominance=1))
    time_final = Sys.time()
    timeMin=as.numeric(difftime(time_final,time_init,units="min"))
    parameters = data.frame(parents=parents,
                            ploidy=ploidy,
                            workingfolder=workingfolder,
                            QTLh2=QTLh2,
                            timeMin=timeMin,
                            sampled_qtlposition=sampleQTL,
                            sampled_qtlname=sampleQTLname)

    output <- list(parameters=parameters,
                   scan=ans,
                   params=list(par1,par2),
                   qtlfitDown=qtlfitDown,
                   qtlfitUp=qtlfitUp,
                   pheno=pheno,
                   SimQTLAlleles=qtlAllelesAll,
                   SimQTLEffects=qtlEffectsAll)
    cat(paste0("Done... it took ",round(timeMin,1)," minutes. \n"))
    setwd(init.wd)
    return(output)
  }else{
    time_final = Sys.time()
    timeMin=as.numeric(difftime(time_final,time_init,units="min"))
    cat(paste0("Done... it took ",round(timeMin,1)," minutes. \n"))
    setwd(init.wd)
    return()
  }
}
