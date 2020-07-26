#'
#' @title Nucleotide Substitution Models for Multiple Alignments
#'
#' @description The function is an extension for multiple alignments of the modelTest() function from phangorn. For more information help("modelTest")
#'
#' @param genes A list objecti with multiple alignments (DNAbin class)
#'
#' @param model The models to evaluate. By default JC, F81, HKY, K80, SYM and GTR. Or use "all". 
#'
#' @param G Boolean, Test Gamma. By default TRUE
#'
#' @param I Boolean, Test invariant sites. By default TRUE
#'
#' @return A list object of length of three. 1. summary.results: A data frame with the selected models for each alignment and each metric (AIC, AICw, AICc, AICcw and BIC). 2. results: Detailed information for each alignment (list object), the normal output of modelTest. 3. warnings: When a alignment have two o more models selected, the warning shows the models selected for each alignment and each metric.  
#'
#' @examples
#'
#' library(alignTools)
#'
#' data(Procyonidae)
#' head(PRocyonidae, 5L)
#'
#' genes <- multiGenBank(Procyonidae, TRUE, FALSE)
#'
#' alignment <- multiMuscle(genes, FALSE)
#'
#' models <- multiModelTest(genes=alignment)
#'
#' models$results.summary
#' str(models$results)
#' models$results$cytb
#' models$warnigs
#'
#' @seealso \code{\link{multiMuscle}}
#'
#' @keywords ModelTest nucleotide substitution
#'
#' @author Omar Daniel Leon-Alvarado <leon.alvarado@@gmail.com>
#'
#' @export

multiModelTest <- function(genes, tree=NULL, model=c("JC","F81","HKY","K80","SYM","GTR"), G=TRUE, I=TRUE, FREQ=FALSE,k=4,
control=pml.control(epsilon=1e-08,maxit=10,trace=1), multicore=FALSE,mc.cores=NULL){

  models <- matrix(NA, nrow=5, ncol=length(genes))
  models.list <- list()
  warnings <- c()

  for(i in 1:length(genes)){
  
    mod.t <- phangorn::modelTest(phangorn::as.phyDat(genes[[i]]), tree=tree, model=model, G=G, I=I, FREQ=FREQ, k=k, control=control, multicore=multicore,mc.cores=mc.cores)

    a <- mod.t$Model[mod.t$AIC==min(mod.t$AIC)]
    if(length(a)>1){
      warnings <- c(warnings,paste0("AIC ",names(genes)[i]," with: ",paste(a,collapse=" ")," models"))
      a <- "*"
    }
    b <- mod.t$Model[mod.t$AICw==min(mod.t$AICw)]
    if(length(b)>1){
      warnings <- c(warnings,paste0("AICw ",names(genes)[i]," with: ",paste(b,collapse=" ")," models"))
      b <- "*"
    }
    c <- mod.t$Model[mod.t$AICc==min(mod.t$AICc)]
    if(length(c)>1){
      warnings <- c(warnings,paste0("AICc ",names(genes)[i]," with: ",paste(c,collapse=" ")," models"))
      c <- "*"
    }
    d <- mod.t$Model[mod.t$AICcw==min(mod.t$AICcw)]
    if(length(d)>1){
      warnings <- c(warnings,paste0("AICcw ",names(genes)[i]," with: ",paste(d,collapse=" ")," models"))
      d <- "*"
    }
    e <- mod.t$Model[mod.t$BIC==min(mod.t$BIC)]
    if(length(e)>1){
      warnings <- c(warnings,paste0("BIC ",names(genes)[i]," with: ",paste(e,collapse=" ")," models"))
      e <- "*"
    }

    models[1,i] <-a
    models[2,i] <-b
    models[3,i] <-c
    models[4,i] <-d
    models[5,i] <-e

    models.list[[i]] <- mod.t

  }

  if(is.null(warnings)){
    warnings <- paste0("No warnings recorded")
  }

  colnames(models) <- names(genes)
  rownames(models) <- c("AIC","AICw","AICc","AICcw","BIC")
  models <- as.data.frame(models)

  names(models.list) <- names(genes)

  warnings <- as.data.frame(warnings)

  output <- list(results.summary=models,results=models.list,warnings=warnings)
  
  print("--------------------------")
  print(models)
  print("If there is a '*', please see the warnings")

 return(output)

}
