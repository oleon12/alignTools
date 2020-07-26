#' @title Write DNA data in TNT format
#'
#' @description Write DNA data in TNT format for one or multiple genes
#'
#' @param genes A list object of multiple alignments (DNAbin class) or a single DNAbin class
#'
#' @param file name of the file. When a list object is used, the file param will be appended to the respective name of each gene.  
#'
#' @examples
#'
#' library(alignTools)
#'
#' data(Procyonidae)
#'
#' genes <- multiGenBank(Procyonidae, TRUE, FALSE)
#' alignment <- multiMuscle(genes, FALSE)
#'
#' #See the  output
#' write.tnt(alignment, file="Procyonidae")
#'
#' write.tnt(alignment[[1]])
#'
#' @seealso \code{\link{read.multiGenes}}
#'
#' @keywords DNA TNT write
#'
#' @author Omar Daniel Leon-Alvarado <leon.alvarado12@@gmail.com>
#'
#' @export

write.tnt <- function(genes, file= NULL, format="dna"){
  
  if(class(genes)!="list"){
    if(length(file)==0){
      filename <- paste0("gene1",".tnt")
    }else{
      filename <- paste0(file,".tnt")
    }
    if(format=="dna"){
      tnt.dna(init=genes, filename=filename)
    }   
  }else{
    for(i in 1:length(genes)){
  
      if(length(file)==0){
        filename <- paste0(names(genes)[i],".tnt")
      }else{
        filename <- paste0(file,names(genes)[i],".tnt")
      }
      if(format=="dna"){
        tnt.dna(init=genes[[i]], filename=filename)
      }
  
    }
  } 
}


