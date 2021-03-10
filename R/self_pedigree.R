#' Build Self Pedigree
#'
#' Build selfs pedigree based on number of parents and population size
#'
#' @param parents number of parents.
#' @param popsize population size.
#' @param padsize string length for the individual name number
#' @param nextinteger if TRUE all the subpopulations will have same size, total population size will changed accordingly to the next possible integer, if FALSE sub populations will be unbalenced
#'
#' @return half-diallel pedigree.
#'
#' @examples
#' ped <- self_pedigree(parents=5,popsize=1000)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export

self_pedigree <- function(parents=NULL,popsize=NULL,padsize=4,nextinteger=FALSE){
  if(parents>1){
    pedigree = diallel_pedigree(parents=parents,subpopsize=popsize/parents,selfs=parents,padsize=padsize)
    pedigree = pedigree[c(1:(popsize+parents)),]
  }else{
    pedigree = diallel_pedigree(parents=(parents+1),subpopsize=popsize,selfs=(parents+1),padsize=padsize)
    pedigree = pedigree[c(1,(3:(popsize+2))),]
  }
  return(pedigree)
}
