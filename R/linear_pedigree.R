#' Build Linear Pedigree
#'
#' Build linear pedigree based on number of parents, number of selfs, and population size
#'
#' @param parents number of parents.
#' @param popsize population size.
#' @param selfs number of selfed populations
#' @param padsize string length for the individual name number
#'
#' @return linear pedigree.
#'
#' @examples
#' ped <- linear_pedigree(parents=4,popsize=100,selfs=1)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export
#'

linear_pedigree = function(parents=4,popsize=100,selfs=1,padsize=4){
  pedigree = round_pedigree(parents=parents,popsize=popsize*parents,selfs=selfs,padsize=padsize)
  pedigree = pedigree[-which(substr(pedigree[,1],1,3)==paste0(LETTERS[parents],"xA")),]
  pedigree$number = substr(pedigree$Name,4,4+padsize)
  pedigree = pedigree[order(pedigree$number),]
  pedigree = pedigree[c(1:(popsize+parents)),-4]
  pedigree = pedigree[order(pedigree$Name),]
  pedigree = pedigree[c(match(LETTERS[1:parents],pedigree$Name),
                  (1:nrow(pedigree))[-match(LETTERS[1:parents],pedigree$Name)]),]
  rownames(pedigree) = NULL
  return(pedigree)
}
