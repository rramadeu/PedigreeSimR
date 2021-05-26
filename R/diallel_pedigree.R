#' Build Diallel Pedigree
#'
#' Build half-diallel pedigree based on number of parents, number of selfs, and population size
#'
#' @param parents number of parents.
#' @param popsize population size.
#' @param subpopsize subpopulation size.
#' @param selfs number of selfed populations
#' @param padsize string length for the individual name number
#' @param nextinteger if TRUE all the subpopulations will have same size, total population size will changed accordingly to the next possible integer, if FALSE sub populations will be unbalenced
#'
#' @return half-diallel pedigree.
#'
#' @examples
#' ped <- diallel_pedigree(parents=5,popsize=1000,selfs=2)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export

diallel_pedigree <- function(parents=NULL,popsize=NULL,subpopsize=NULL,selfs=0,padsize=4,nextinteger=FALSE){

  ## Checks
  if(selfs>parents)
    stop(deparse("selfs should be equal or less than parent"))

  if(!is.null(popsize) && !is.null(subpopsize))
    stop(deparse("popsize or subpopsize should be NULL, just indicate one"))

  if(parents>702)
    stop(deparse("parents should be lower than 702"))
  tmp <- expand.grid(LETTERS, LETTERS) #expanding for more than 26 parents
  tmp <- tmp[order(tmp$Var1,tmp$Var2),] #expanding for more than 26 parents
  LETTERS2 <- c(LETTERS, do.call('paste0',tmp)) #expanding for more than 26 parents
  
  # Build unknown pedigree
  pedigree = data.frame(Name = LETTERS2[1:parents],
                        Parent1 = "NA",
                        Parent2 = "NA")

  subpops = ncol(combn(parents,2))+selfs

  # Defining subpopsize and popsize
  if(is.null(subpopsize)){
    subpopsize = popsize/subpops
  }else{
    popsize=subpopsize*subpops
  }

  origpopsize=popsize

  # Rounding subpopsize to the next possible integer
  while(subpopsize %% floor(subpopsize) > 0){
    popsize=popsize+1
    subpopsize = (popsize)/subpops
  }

  if(padsize<ceiling(log(subpopsize,10))+1)
    padsize=ceiling(log(subpopsize,10))+1

  f1code = stringr::str_pad(1:subpopsize,
                            width = padsize,
                            pad = 0)

  totalsubpop = popsize/subpopsize
  if(nextinteger==FALSE){
    subpopsize = rep(subpopsize,totalsubpop)-c(rep(0,totalsubpop-(popsize-origpopsize)),rep(1,popsize-origpopsize))
  }else{
    subpopsize = rep(subpopsize,totalsubpop)
  }

  # Build selfs pedigree
  if( selfs>0){
    f1codeselfs=NULL
    for(i in 1:selfs)
      f1codeselfs = c(f1codeselfs,f1code[1:subpopsize[i]])
    selfsped = data.frame(Name = paste0(rep(LETTERS2[1:selfs],times=subpopsize[1:selfs]),
                                        "x",
                                        rep(LETTERS2[1:selfs],times=subpopsize[1:selfs]),
                                        f1codeselfs),
                          Parent1 = rep(LETTERS2[1:selfs],times=subpopsize[1:selfs]),
                          Parent2 = rep(LETTERS2[1:selfs],times=subpopsize[1:selfs]))
  }else{
    selfsped=NULL
  }

  if(selfs>0)
    subpopsize=subpopsize[-c(1:selfs)]
  f1codediallel = NULL
  for(i in 1:ncol(combn(parents,2)))
    f1codediallel = c(f1codediallel,f1code[1:subpopsize[i]])
  Parent1 = LETTERS2[rep(combn(parents,2)[1,],times=subpopsize)]
  Parent2 = LETTERS2[rep(combn(parents,2)[2,],times=subpopsize)]
  Name = paste0(Parent1,"x",Parent2,f1codediallel)
  diallelped = data.frame(Name=Name,Parent1=Parent1,Parent2=Parent2)

  return(pedigree = rbind(pedigree,selfsped,diallelped))
}



