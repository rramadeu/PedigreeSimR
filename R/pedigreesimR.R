#' PedigreeSim Management
#'
#' Wrap-up function to run PedigreeSim software
#'
#' @param map vector of the marker position.
#' @param haplotypes matrix with haplotype information to be draw
#' @param pedigree pedigree data frame
#' @param centromere numeric with its position
#' @param prefPairing numeric
#' @param quadrivalents numeric
#' @param ploidy integer even number
#' @param output string with the prefix of the output fule
#' @param mapfunction "HALDANE" or "KOSAMBI"
#' @param chromosome string
#' @param seed integer to be used to sample haplotypes
#' @param allownochiasmata numeric
#' @param naturalpairing numeric
#' @param parallelquadrivalents numeric
#' @param pairedcentromeres numeric
#' @param mapwidthpad numeric length of marker code
#'
#' @return nothing?
#'
#' @examples
#' map = 1:100
#' haplotypes = fake_haplo(100,100,seed=1234)
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
                         ploidy=2,
                         workingfolder="PedigreeSimR_files",
                         mapfunction="HALDANE",
                         chromosome="A",
                         seed=NULL,
                         allownochiasmata=1,
                         naturalpairing=1,
                         parallelquadrivalents=0,
                         pairedcentromeres=0,
                         mapwidthpad=4){

  ## Creating map file
  mapdf = data.frame(marker=paste0(chromosome,"_",str_pad(map,width = mapwidthpad,side = "left",pad = "0")),
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
  haplotypes = haplotypes[,sampledhaplos]
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
                                      paste0(workingfolder,"/input.chrom"),
                                      paste0(workingfolder,"/input.ped"),
                                      paste0(workingfolder,"/input.map"),
                                      paste0(workingfolder,"/input.founder"),
                                      paste0(workingfolder,"/out."),
                                      allownochiasmata,
                                      naturalpairing,
                                      parallelquadrivalents,
                                      pairedcentromeres)
                             )

  ## Writing files
  if(is.na(match(workingfolder,list.files())))
    system(paste0("mkdir ",workingfolder))

  write.table(chrdf,file=paste0(workingfolder,"/input.chrom"),row.names = FALSE,quote = FALSE)
  write.table(pedigree,file=paste0(workingfolder,"/input.ped"),row.names = FALSE,quote = FALSE)
  write.table(mapdf,file=paste0(workingfolder,"/input.map"),row.names = FALSE,quote = FALSE)
  write.table(founderdf,file=paste0(workingfolder,"/input.founder"),row.names = FALSE,quote = FALSE)
  write.table(parameterdf,file=paste0(workingfolder,"/input.par"),row.names = FALSE, col.names=FALSE,quote = FALSE)

  system(
    paste0("java -jar ",system.file(package = "PedigreeSimR"),"/PedigreeSim2/PedigreeSim.jar ",workingfolder,"/input.par")
  )

}
