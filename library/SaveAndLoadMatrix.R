
getColnames <- function(path)
{
  #get filename
  colsFile <- paste(path,"barcodes.tsv",sep = "/")
  #read file
  cols<-read.csv(colsFile,header = FALSE)
  #select column
  cols<-c(cols$V1)
  #return
  cols
}

getRowNames <- function(path)
{
  #get filename
  rowFile  <- paste(path,"features.tsv",sep = "/")
  #read file
  rows <- read.csv(rowFile,sep = "\t",header = FALSE)
  #select column
  rows<-c(rows$V2)
  #return
  rows
}

loadMatrix <- function(path)
{
  rows <- getRowNames(path)
  cols <- getColnames(path)
  
  library(Matrix)
  
  matrixFile<-paste(path,"matrix.mtx",sep="/")
  
  #read matrix file
  matrix <- readMM(matrixFile)
  
  #set rows and cols names
  rownames(matrix)<-rows
  colnames(matrix)<-cols
  
  #return
  matrix
}




getGenePositions <- function(matrix,gene)
{
  
  cells_indexes <- read.csv("res/cells_indexes.csv")
  #order and sort positions
  cols <- colnames(matrix)

  positions <- cells_indexes[match(cols,cells_indexes$barcode),]
  
  #initialize variables
  x<-c(NULL)
  y<-c(NULL)
  val<-c(NULL)
  #extract gene expression vector
  gene_expr<-matrix[gene,]
  
  #for each row take gene expression
  for (i in 1:length(gene_expr)) {
    if(gene_expr[i] > 0)
    {
      val<-append(val,gene_expr[i])
      x<-append(x,positions$array_col[i])
      y<-append(y,positions$array_row[i])
    }
  }
  
  as.data.frame(cbind(x,y,val))
}

