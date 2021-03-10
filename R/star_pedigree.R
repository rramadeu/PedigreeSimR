#' Build Star Pedigree
#'
#' Build star pedigree based on number of parents, number of selfs, and population size
#'
#' @param parents number of parents.
#' @param popsize population size.
#' @param selfs number of selfed populations
#' @param padsize string length for the individual name number
#'
#' @return star pedigree.
#'
#' @examples
#' ped <- star_pedigree(parents=4,popsize=100,selfs=1)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export
#'
star_pedigree = function(parents=4,popsize=100,selfs=1,padsize=4){
  pedigree = diallel_pedigree(parents=parents,popsize=popsize*parents,selfs=selfs,padsize=padsize)
  pedigree = pedigree[unique(c(which(as.character(pedigree[,2])==as.character(pedigree[,3])),
                               which(pedigree[,2] %in% c("A")))),]
  pedigree$number = substr(pedigree$Name,4,4+padsize)
  pedigree = pedigree[order(pedigree$number),]
  pedigree = pedigree[c(1:(popsize+parents)),-4]
  pedigree = pedigree[order(pedigree$Name),]
  pedigree = pedigree[c(match(LETTERS[1:parents],pedigree$Name),
                      (1:nrow(pedigree))[-match(LETTERS[1:parents],pedigree$Name)]),]
  rownames(pedigree) = NULL
  return(pedigree)
}
