library(Matrix)
library(ggplot2)
library(tidymodels)
library(xgboost)
library(edgeR)
set.seed(12345)

#HUMAN BREAST CANCER FROM 10XGENOMICS VISIUM SPATIAL PROTOCOLS, Luminal B

spatial_expression <- as.matrix(readMM("data/spatial/matrix.mtx"))

#NORMALIZATION
spatial_expression<- cpm(calcNormFactors(DGEList(counts=spatial_expression),normalized.lib.sizes = T))
spatial_expression<-log10(spatial_expression+1)
spatial_expression <- t(spatial_expression)


features <- read.csv("data/spatial/features.tsv", sep="\t",header = F)
barcodes <- read.csv("data/spatial/barcodes.tsv",sep="\t",header =F )
indexes <- read.csv("data/spatial/cells_indexes.csv") 

indexes <- indexes[match(barcodes$V1,indexes$barcode),]

DATA <- data.frame(spatial_expression)
DATA <- cbind(indexes,DATA)

#removing barcodes with few reads
percentage <- (ncol(DATA)/100)*(95)
DATA = DATA[rowSums(DATA[,4:36604]==0)<=percentage,]
colnames(DATA)[4:36604] <- features$V1

barcodes <- DATA[1:3]
pam50 <- read.csv("resources/pam50.csv")
pam50_reads <- DATA[pam50$ens]

remove(DATA,pam50,percentage,features,indexes,spatial_expression)


# let's predict

model <- readRDS("result/model_xgboost.rds")

res <- predict(model,pam50_reads)

res <- cbind(barcodes,res)

ggplot(data=res,aes(x=array_col,y=array_row, colour=.pred_class))+geom_point()

# mainly LUMINAL B !

