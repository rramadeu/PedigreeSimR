#' Build fake haplotypes
#'
#' Build fake haplotypes based on given number of markers and number of haplotypes
#'
#' @param n number of haplotypes (integer).
#' @param m number of markers (integer).
#' @param seed seed to be used in the sampling
#'
#' @return Matrix with the Relationship between the individuals.
#'
#' @examples
#' haps <- fake_haplo(n=10,m=100,seed=12345)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export

fake_haplo <- function(n=NULL,m=NULL,seed=NULL){
  if(!is.null(seed))(set.seed(seed))
  haplo <- sample(c(0,1),n*m,TRUE)
  haplo <- matrix(haplo,ncol=n)
  return(haplo)
}


