#'
#' @title TNT format syntax (DNA)
#'
#' @param init DNAbin object
#'
#' @param filename output file name
#'
#' @keywords DNA TNT write
#'
#' @author Omar Daniel Leon-Alvarado <leon.alvarado12@@gmail.com>
#' 
#' @export

tnt.dna <- function(init, filename){

  nb.l <- dim(init)[2]
  sp.l <- dim(init)[1]

  seq <- as.alignment(init)

  cat("xread", file = filename)
  cat("\n", file = filename, append = T)
  cat(c(nb.l,sp.l), file = filename, append = T)
  cat("\n", file = filename, append = T)
  cat("&[ dna ]", file = filename, append = T)
  cat("\n", file = filename, append = T)
  
  for(x in 1:length(seq$nam)){
    cat(paste(seq$nam[x],seq$seq[x], sep=" "),file=filename, append=T)
    cat("\n", file = filename, append = T)
  }

  cat(";", file = filename, append = T)
  cat("\n", file = filename, append = T)
  cat("proc/;", file = filename, append = T)


}
