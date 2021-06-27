


gene_id_to_symbol <- function (id){
  
  dictionary <- read.csv("res/dictionary.csv")
  dictionary[dictionary$ID==id,2]
}

gene_symbol_to_id <- function (name){
  dictionary[dictionary$symbol==name,1]
}

removeUnexpressed <- function(matrix,percentage){

  percentage <- (ncol(matrix)/100)*(100-percentage)
  matrix[rowSums(matrix==0)<=percentage,]
}


logMatrix <- function(matrix){
  log10(matrix+1.0)
}

filter_gene_count_matrix <- function (matrix,filter) {

  matrix[names(matrix)[names(matrix) %in% filter]]
}