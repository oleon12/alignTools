#'
#' @title Concatenate Multiple Alignments
#'
#' @description This function concatenate multiple alignments
#'
#' @param genes A list object with multiple alignments (DNAbin class)
#'
#' @param missing When a species or individual have no information for one or more genes, the absent gene(s) should be filled with "?" or "-". By default "?".
#'
#' @param write.dna Boolean, write the concatenate alignment in a file with predeterminate format. By deafault TRUE.
#'
#' @param write.format The format to write the concatenate alignment (only if write.dna == TRUE). Three options: 1. nexus (default), 2. phylip, 3. fasta or 4. tnt. 
#'
#' @param filename The name of the output file (only if write.dna==TRUE).
#'
#' @return A concatenate alignment of DNAbin class
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
#' alignment <- multiMuscle(x=genes, FALSE)
#'
#' concat <- concatGenes(genes=alignment, missing="?", write.dna=TRUE, write.format="phylip", filename="Concat_eg_output")
#'
#' str(concat)
#'
#' @seealso \code{\link{multiMuscle}} \code{\link{completeGenes}}
#'
#' @keywords concatenate genes alignments
#'
#' @author Omar Daniel Leon-Alvarado <leon.alvarado12@@gmail.com>
#'
#' @includeRmd vignettes/alignTools.rmd
#'
#' @export

concatGenes <- function(Genes, missing =c("?","-"), write.dna=c(TRUE,FALSE), write.format=c("nexus","phylip","fasta","tnt"),filename=NULL){

  #Control1
  if(class(Genes)!="list"){
   stop("Genes must be a list object")
  }

  #Control2
  if(length(missing)==2){

   miss <- "?"

  }
  if(length(missing)==1){

  miss <- missing

  }
  #Control3
  if(!missing%in%c("?","-")){
    stop("Please set an appropriate 'missing' paramenter = ? or -")
  }

  #Control4
  if(length(write.dna)==2){
    write.dna <- TRUE
  }

  pb <- txtProgressBar(min = 0, max = (length(Genes)-1), style = 3)

  for(i in 1:(length(Genes)-1)){

  if(i == 1){
  
    lseq1 <- dim(Genes[[i]])[2]
    lseq2 <- dim(Genes[[i+1]])[2]
  
    seq1 <- ape::as.alignment(Genes[[i]])
    seq2 <- ape::as.alignment(Genes[[i+1]])
  
    spAll <- unique(c(seq1$nam, seq2$nam))
  
    seqs <- c()
  
    for(j in 1:length(spAll)){
    
      pos1 <- which(seq1$nam %in% spAll[j])
      pos2 <- which(seq2$nam %in% spAll[j])
    
      if(length(pos1) != 0){
        seqind1 <- seq1$seq[pos1]
      }else{
        seqind1 <- rep(miss, lseq1)
        seqind1 <- paste(seqind1, collapse = '')
      }
    
      if(length(pos2) != 0){
        seqind2 <- seq2$seq[pos2]
      }else{
        seqind2 <- rep(miss, lseq2)
        seqind2 <- paste(seqind2, collapse = '')
      }
    
      seqTMP <- paste(seqind1, seqind2, sep="")
    
      seqs <- c(seqs, seqTMP)
    
    }
  
    bin1 <- list(nb = length(spAll),
                 seq = seqs,
                 nam = spAll,
                 com = NA)
    class(bin1) <- "alignment"
  
    bin1 <- ape::as.DNAbin(bin1)
  
  }

  if(i != 1){
  
    lseq1 <- dim(bin1)[2]
    lseq2 <- dim(Genes[[i+1]])[2]  

    seq1 <- ape::as.alignment(bin1)
    seq2 <- ape::as.alignment(Genes[[i+1]])
  
    spAll <- unique(c(seq1$nam, seq2$nam))
  
    seqs <- c()
  
    for(j in 1:length(spAll)){
    
      pos1 <- which(seq1$nam %in% spAll[j])
      pos2 <- which(seq2$nam %in% spAll[j])
    
      if(length(pos1) != 0){
        seqind1 <- seq1$seq[pos1]
      }else{
        seqind1 <- rep(miss, lseq1)
        seqind1 <- paste(seqind1, collapse = '')
      }
    
      if(length(pos2) != 0){
        seqind2 <- seq2$seq[pos2]
      }else{
        seqind2 <- rep(miss, lseq2)
        seqind2 <- paste(seqind2, collapse = '')
      }
    
      seqTMP <- paste(seqind1, seqind2, sep="")
    
      seqs <- c(seqs, seqTMP)
    
    }
  
    bin1 <- list(nb = length(spAll),
                 seq = seqs,
                 nam = spAll,
                 com = NA)
    class(bin1) <- "alignment"
  
    bin1 <- ape::as.DNAbin(bin1)  
    }
  
  setTxtProgressBar(pb, i)

  }

  close(pb)

  if(write.dna==TRUE){
    if(length(write.format)==4){
      write.format <- "nexus"
    }
    if(!write.format%in%c("nexus","phylip","fasta","tnt")){
      stop("Please set an appropriate write.format")
    }
    if(is.null(filename)){
      stop("A file name is required")
    }

    if(write.format=="nexus"){
      ape::write.nexus.data(bin1, file=paste(filename,"nex",sep="."), format="dna", interleaved=F)
    }
    if(write.format=="phylip"){
      phangorn::write.phyDat(phangorn::as.phyDat(bin1), format="sequential", file=paste(filename,"phy",sep="."))
    }
    if(write.format=="fasta"){
      ape::write.dna(bin1, file=paste(filename,"fas",sep="."), format="fasta")
    }
    if(write.format=="tnt"){
      write.tnt(bin1, file=filename, format="dna")
    }
  }

print(partitionScheme(genes=Genes))

return(bin1)

}
