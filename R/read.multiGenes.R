#'
#' @title Read Multiple Genes/Alignments files
#'
#' @description The function reads two o more dna files (alignments)
#'
#' @param files Vecotr with the name of the files to read. The files must be in the same working directory.
#'
#' @param format The format of the files. Three options: 1. nexus, 2. phylip or 3. fasta
#'
#' @param names (optional), a vector with same lenght as files. This represent an specific name for each file. If missin, the function will use the files vector as names. 
#'
#' @examples
#' 
#' library(alignTools)
#'
#' input <- c("files1.nex","file2.nex","file3.nex")
#'
#' # In this example, names will be the same as input
#' data <- read.multiGenes(files=input, format="nexus")
#'
#' # In this exmaple, names will be "id1","id2","id3".
#' data <- read.multiGenes(files=input, format="nexus",names=c("id1","id2","id3"))
#'
#' str(data)
#'
#' @seealso \code{\link{multiGenBank}} \code{\link{concatGenes}}
#'
#' @keywords read multiple genes alignments
#'
#' @author Omar Daniel Leon-Alvarado
#'
#' @export

read.multiGenes <- function(files, format=c("nexus","phylip","fasta"), names=NULL){
  #Control1
  if(length(names)!=0){
    if(length(names)!=length(files)){
      stop("Files and names must have the same length")
    }
  }else{
    names <- files
  }

  #Control2
  if(length(format)==3){
    stop("Please select one format: nexus, phylip or fasta")
  }
  
  r.list <- c()
  
  if(!format%in%c("nexus","phylip","fasta")){
    stop("Please set an appropiate format")
  }

  for(i in 1:length(files)){
    if(format=="nexus"){
      r.list[[i]] <- ape::as.DNAbin(ape::read.nexus.data(files[i]))
    }
    if(format=="phylip"){
      r.list[[i]] <- ape::as.DNAbin(phangorn::read.phyDat(files[i]))
    }
    if(format=="fasta"){
      r.list[[i]] <- ape::read.dna(files[i], format="fasta")
    }
  }
  
  names(r.list) <- names
  return(r.list)
  
}
