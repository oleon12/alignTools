#'
#' @title Multiple Genes Alignment with MUSCLE
#'
#' @description This function is an extesion the muscle() function from ape. The x parameter is the only changed here, so, for more information about the other parameters, please see the documentation from ape.
#'
#' @param x A list object with multiple DNAbin objects
#'
#' @param write.dna Boolean, write each alignment as separate files in a predeterminate format. The file names will be the same as the x object names. By deafault TRUE. 
#'
#' @param write.format The format to write each alignment (only if write.dna == TRUE). Four options: 1. nexus (Default), 2. phylip, 3. fasta or 4. tnt. 
#'
#' @return A list object with multiple alignments (DNAbin class). 
#'
#' @examples
#'
#' library(alignTools)
#'
#' data(Procyonidae)
#' head(Procyonidae, 5L)
#'
#' genes <- multiGenBank(Procyonidae, TRUE, FALSE)
#' 
#' alignment <- multiMuscle(x=genes, write.dna=TRUE, write.format="phylip")
#'
#' str(alignment) 
#' alignment$rag2
#'
#' @seealso \code{\link{multiGenBank}}
#'
#' @keywords alignment muscle genes
#'
#' @author Omar Daniel Leon-Alvarado <leon.alvarado12@@gmail.com>
#'
#' @export

multiMuscle <- function(x=NULL,y, guide.tree, quiet=TRUE, original.ordering=TRUE, MoreArgs="", write.dna=c("TRUE","FALSE"), write.format=c("nexus","phylip","fasta","tnt")){

  #Control1
  if(length(write.dna)==2){
    write.dna <- TRUE
  }

  #Control2
  if(is.null(names(x))){
    stop("x need names, set it with names() function")
  }
  #Control3
  if(any(is.na(names(x)))){
   stop("Names of x can not be NA")
  }


  out <- list()

  pb <- txtProgressBar(min = 0, max = length(x), style = 3)

  for(i in 1:length(x)){

    setTxtProgressBar(pb, i)

    if(missing(y)){
      if(missing(guide.tree)){
        alig <- ape::muscle(x=x[[i]], quiet=quiet, original.ordering=TRUE, MoreArgs=MoreArgs) 
      }else{
         alig <- ape::muscle(x=x[[i]], guide.tree=guide.tree, quiet=quiet, original.ordering=original.ordering,MoreArgs=MoreArgs)
       }
    }else{
       alig <- ape::muscle(x=x[[i]],y=y,guide.tree=guide.tree,quiet=quiet,original.ordering=original.ordering,MoreArgs=MoreArgs) 
     }

    out[[i]] <- alig

  }

  close(pb)

  names(out) <- names(x)

  if(write.dna==TRUE){
    if(length(write.format)==4){
      write.format <- "nexus"
    }
    if(!write.format%in%c("nexus","phylip","fasta","tnt")){
      stop("Please set an appropriate write.format")
    }

    for(i in 1:length(out)){
      if(write.format=="nexus"){
        ape::write.nexus.data(out[[i]],file=paste(names(out)[i],"nex",sep="."), format="dna", interleaved=F)
      }
      if(write.format=="phylip"){
        phangorn::write.phyDat(phangorn::as.phyDat(out[[i]]), format="sequential", file=paste(names(out)[i],"phy",sep="."))
      }
      if(write.format=="fasta"){
        ape::write.dna(out[[i]], file=paste(names(out)[i],"fas",sep="."), format="fasta")
      }
      if(write.format=="tnt"){
        write.tnt(out[[i]], file=names(out)[i], format="dna")
      }
    }
  }

  return(out)

}
