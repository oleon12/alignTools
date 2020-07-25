#'
#' @title Complete Alignments
#'
#' @description Some softwares (e.g. RevBayes) for partitioned analysis required each gene in a separate file, but all files must have the same number of species. This function achive that task. 
#'
#' @param genes A list object with multiple alignments (DNAbin class)
#'
#' @param missing When a species or individual have no information for a gene, the absent gene should be filled with "?" or "-". by default "?".
#'
#' @param write.dna Boolean, write on separate files each gene with the complete list of species or individuals. By default TRUE.
#'
#' @param write.format The format to write each alignment (only if write.dna==TRUE). Four options: 1. nexus (default), 2. phylip, 3. fasta or 4. tnt.
#'
#' @return A list object with multiple alignments (DNAbin class). All with the same number of species or individuals.
#'
#' @examples
#'
#' library(alignTools)
#'
#' data(Procyonidae)
#' head(Procyonidae, 5L)
#'
#' genes <- multiGenBank(Procyonidae,TRUE,FALSE)
#'
#' alignment <- multiMuscle(genes, FALSE)
#'
#' complete <- completeGenes(genes=alignment, missing="?", write.dna=TRUE, write.format="fasta")
#'
#' str(complete)
#'
#' @seealso \code{\link{concatGenes}} \code{\link{read.multiGenes}}
#'
#' @keywords complete genes alignments
#'
#' @author Omar Daniel Leon-Alvarado <leon.alvarado12@@gmail.com>
#'
#' @includeRmd vignettes/alignTools.rmd
#'
#' @export

completeGenes <- function(genes, missing=c("?","-"),write.dna=c(TRUE,FALSE), write.format=c("nexus","phylip","fasta","tnt")){
  
  #Contro1
  if(length(missing)==2){
    missing <- "?"
  }
  #Control2
  if(length(write.dna)==2){
    write.dna <- TRUE
  }
  
  all <- as.alignment(concatGenes(genes, missing = missing,write.dna = FALSE))$nam
  out.list <- list()
  
  pb <- txtProgressBar(min = 0, max = length(genes), style = 3)

  for(i in 1:length(genes)){
    in.seq <- as.alignment(genes[[i]])
    lseq <- dim(genes[[i]])[2]
    
    miss <- all[-which(all%in%in.seq$nam)]
    
    if(length(miss)!=0){
      for(j in 1:length(miss)){
        in.seq$seq[in.seq$nb+j] <- paste(rep(missing, lseq), collapse = "")
        in.seq$nam[in.seq$nb+j] <- miss[j]
      }
      in.seq$nb <- length(in.seq$nam)
      out.list[[i]] <- as.DNAbin(in.seq)
    }else{
      out.list[[i]] <- genes[[i]]
    }
    
    if(write.dna==TRUE){
      if(length(write.format)==4){
        write.format  <- "nexus"
      }
      if(!write.format%in%c("nexus","phylip","fasta","tnt")){
        stop("Please set an appropiate write.format")
      }
      if(write.format=="nexus"){
        ape::write.nexus.data(out.list[[i]], file=paste0(names(genes)[i],".nex"), format="dna", interleaved=F)
      }
      if(write.format=="phylip"){
        phangorn::write.phyDat(phangorn::as.phyDat(out.list[[i]]), format="sequential", file=paste0(names(genes)[i], ".phy"))
      }
      if(write.format=="fasta"){
        ape::write.dna(out.list[[i]], file=paste0(names(genes)[i],".fas"), format="fasta")
      }
      if(write.format=="tnt"){
        write.tnt(out.list[[i]], file=names(genes)[i], format="dna")
      }
    }
 
  setTxtProgressBar(pb,i)

  }

  close(pb)
  
  names(out.list) <- names(genes)
  return(out.list)
  
}
