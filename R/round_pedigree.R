#' Build Round-Robin Pedigree
#'
#' Build single round-robin pedigree based on number of parents and population size
#'
#' @param parents number of parents.
#' @param popsize population size.
#' @param subpopsize subpopulation size.
#' @param padize string length for the individual name number
#' @param nextinteger if TRUE all the subpopulations will have same size, total population size will changed accordingly to the next possible integer, if FALSE sub populations will be unbalenced
#'
#' @return single round-robin pedigree.
#'
#' @examples
#' haps <- fake_haplo(n=10,map=100,seed=12345)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export


round_pedigree <- function(parents=NULL,popsize=NULL,subpopsize=NULL,padsize=4,nextinteger=FALSE){

  if(!is.null(popsize) && !is.null(subpopsize))
    stop(deparse("popsize or subpopsize should be NULL, just indicate one"))

  pedigree = data.frame(Name = LETTERS[1:parents],
                        Parent1 = "NA",
                        Parent2 = "NA")

  subpops = parents

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


  f1coderobin = NULL
  for(i in 1:subpops)
    f1coderobin = c(f1coderobin,f1code[1:subpopsize[i]])
  Parent1 = LETTERS[rep(1:subpops,times=subpopsize)]
  Parent2 = LETTERS[rep(c(2:subpops,1),times=subpopsize)]
  Name = paste0(Parent1,"x",Parent2,f1coderobin)
  robinped = data.frame(Name=Name,Parent1=Parent1,Parent2=Parent2)

  return(pedigree = rbind(pedigree,robinped))
}
