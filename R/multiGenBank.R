#' 
#' @title Get Multiple Sequences from GenBank
#'
#' @description The function download the sequences for multiple genes. The sequences will be stored in a list object and/or will be saved in different files in fasta format.
#'
#' @param x is the matrix or data frame with the access numbers of each gene for each species or individual. This matrix must have the species or indivual names in the first column. See example.
#'
#' @param with.na Boolean, if there is any NA in the matrix, must be set as TRUE (default). When you have no data or access number in a gene for one or more species, you must codify this as NA, otherwise, the function will not work. 
#'
#' @param write.dna Boolean, write each gene as separate files in fasta format. The file names will be the same as the column names. By deafault TRUE. 
#'
#' @return The function returns a list object with the sequences for each gene (DNAbin class). Even setting write.dna=TRUE, the function will returns the list object.  
#'
#' @examples
#'
#' library(alignTools)
#'
#' data(Procyonidae)
#' head(Procyonidae, 5L)
#'
#' genes <- multiGenBank(x=Procyonidae, with.na=TRUE, write.dna=FALSE)
#' str(genes)
#' genes$cytb
#'
#' @seealso \code{\link{multiMuscle}} \code{\link{read.multiGenes}}
#'
#' @keywords GenBank sequences download
#'
#' @author Omar Daniel Leon-Alvarado <leon.alvarado12@@gmail.com>
#'
#' @includeRmd vignettes/alignTools.rmd 
#'
#' @export

multiGenBank <- function(x, with.na=c("TRUE","FALSE"), write.dna=c("TRUE","FALSE")){

  #Control1
  if(length(with.na)==2){
    with.na <- TRUE
  }
  #Control2
  if(length(write.dna)==2){
    write.dna <- TRUE
  }

  ####
  genes <- x[ ,2:length(colnames(x))]
  genes.n <- colnames(x)[2:length(colnames(x))]
  gene.list <- list()

  ####
  # Only one gene
  ####

  if(length(colnames(genes))==0){
    gen <- genes
    names(gen) <- x[,1]
    gen.nm <- names(gen)

    if(with.na==TRUE){
      gen.nm <- names(na.omit(gen))
      gen <- na.omit(gen)
    }

    seq <- ape::read.GenBank(gen)
    names(seq) <- gen.nm
  
    gene.list[[1]] <- seq
    names(gene.list) <- genes.n

    if(write.dna==TRUE){
      file.id <- paste(names(genes.list)[[1]], ".fas", sep="")
      ape::write.dna(gene.list[[1]], file=file.id, format="fasta")
    }
    
    return(gene.list)
  }


  ######
  #Only for multiple genes
  ######

  if(with.na==TRUE){
    pb <- txtProgressBar(min=0, max=length(colnames(genes)),style=3)

    for(i in 1:length(colnames(genes))){
      gen <- genes[,i]
      names(gen) <- x[,1]
      
      gen <- na.omit(gen)

      seq <- ape::read.GenBank(gen)
      names(seq) <- names(gen)

      gene.list[[i]] <- seq

      setTxtProgressBar(pb,i)

    }

    close(pb)
    names(gene.list) <- genes.n
 
  }

  if(with.na==FALSE){
    pb <- txtProgressBar(min=0, max=length(colnames(genes)),style=3)

    for(i in 1:length(colnames(genes))){
      gen <- genes[,i]
      #Control3
      if(length(which(is.na(gen)))!=0){
        stop("The data have NA, please set with.na TRUE")
      }
      names(gen) <- x[,1]
      seq <- ape::read.GenBank(gen)
      names(seq) <- names(gen)

      gene.list[[i]] <- seq

      setTxtProgressBar(pb,i)
    
    }

    close(pb)
    names(gene.list) <- genes.n
  }

  if(write.dna == TRUE){
    for(i in 1:length(gene.list)){
      file.id <- paste(names(gene.list)[i], ".fas", sep="")
      ape::write.dna(gene.list[[i]], file=file.id,format="fasta")
    }
  }

  return(gene.list)
}
