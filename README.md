# alignTools: easy download and manage GenBank data, and alignments for phylogenetics </br>
**Author:** Omar Daniel Leon-Alvarado, leon.alvarado12@gmail.com
***
## **Installation**
```{r}
library(devtolls)
install_github("oleon12/alignTools")
```
***
## **Data**
To download multiple genes from GenBank, the input file we use is a matrix. This matrix should have in the first column the species names or identifiers, and the others columns will be the genes and their accession numbers. So, this is the "hard" step, you need to find each accession number. When a species have not information for a gene, you fill it as NA. See the Procyonidae data.  
```{r}
#Procyonidae data
head(Procyonidae, 5L)
```

