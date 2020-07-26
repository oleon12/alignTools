#'
#' @title Partition Scheme for a Concatenate Alignment
#'
#' @description Shows a partition scheme for a concatenate alignment of N genes
#'
#' @param genes A list objecto with multiple alginments (DNAbin class)
#'
#' @examples
#'
#' library(alignTools)
#'
#' data(Procyonidae)
#' head(Procyonidae, 5L)
#'
#' genes <- multiGenBank(Procyonidae,TRUE,FALSE)
#' alignment <- multiMuscle(genes, write.dna=FALSE)
#' 
#' concat <- concatGenes(alignment, "?", FALSE)
#'
#' partitionScheme(alignment)
#'
#' @seealso \code{\link{concatGenes}} \code{\link{completeGenes}}
#'
#' @keywords partition scheme genes alignments
#'
#' @author Omar Daniel Leon-Alvarado 
#'
#' @export

partitionScheme <- function(genes){

  le <- unlist(lapply(genes, function(x){dim(x)[2]}))
  le.p <- c()

  for(i in 1:length(le)){
    if(i==1){
      l2 <- c(0,le[-length(le)])
      le.p <- le+l2
    }else{
      l2 <- c(0,l2[-length(l2)])
      le.p <- le.p+l2
    }
  }

  pos1 <- le.p+1
  pos1 <- c(1, pos1[-length(pos1)])

  out <- paste(pos1, le.p, sep="-")
  names(out) <- names(genes)

  print("Partition Scheme")
  return(out)
}
