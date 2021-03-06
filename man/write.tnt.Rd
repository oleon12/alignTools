% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write.tnt.R
\name{write.tnt}
\alias{write.tnt}
\title{Write DNA data in TNT format}
\usage{
write.tnt(genes, file = NULL, format = "dna")
}
\arguments{
\item{genes}{A list object of multiple alignments (DNAbin class) or a single DNAbin class}

\item{file}{name of the file. When a list object is used, the file param will be appended to the respective name of each gene.}
}
\description{
Write DNA data in TNT format for one or multiple genes
}
\examples{

library(alignTools)

data(Procyonidae)

genes <- multiGenBank(Procyonidae, TRUE, FALSE)
alignment <- multiMuscle(genes, FALSE)

#See the  output
write.tnt(alignment, file="Procyonidae")

write.tnt(alignment[[1]])

}
\seealso{
\code{\link{read.multiGenes}}
}
\author{
Omar Daniel Leon-Alvarado \href{mailto:leon.alvarado12@gmail.com}{leon.alvarado12@gmail.com}
}
\keyword{DNA}
\keyword{TNT}
\keyword{write}
