source("library/utils.R")
library(ggplot2)
library(tidymodels)
library(xgboost)
set.seed(12345)


########################DATA LOADING
dataset <- NA

if(!file.exists("data/cancer-subtype-seqdata.csv")){ # load data from tcga
  
  source("library/tcga-brca.R")
  
  experiment <- load_brca_experiment_data()
  clinic_data <- as.data.frame(get_brca_clinic_data(experiment))
  gene_expr <- get_brca_genes_count(experiment)
  
  Subtype <-clinic_data$paper_BRCA_Subtype_PAM50
  
  
  ds<- cbind(Subtype,gene_expr)
  ds <- na.omit(ds)
  ds <- ds[ds$Subtype!="Normal",]
  write.csv(ds,"data/cancer-subtype-seqdata.csv")
}

############################### PREPROCESSING

if(!file.exists("data/brca_tcga_pam50.csv")) {  # normalize data and get gene  signature
  
  dataset <- read.csv("output/cancer-subtype-seqdata.csv",row.names = 1)
  Subtype <- dataset$Subtype
  
  matrix <- dataset[-1]
  
  library(edgeR)
  source("library/pam50.R")
  source("library/geneOperation.R")
  
  matrix <- t(matrix)
  matrix<- cpm(calcNormFactors(DGEList(counts=matrix),normalized.lib.sizes = T))
  matrix<-log10(matrix+1)
  matrix <- t(matrix)
  
  
  dataset <- cbind(Subtype,matrix[get_pam50_ids()])
  
  write.csv(dataset,"output/brca_tcga_pam50.csv")
  
}else{
  
  dataset <- read.csv("data/brca_tcga_pam50.csv",row.names = 1, stringsAsFactors = T)
  
}


class_numerosity = data.frame(table(dataset$Subtype))
ggplot(data=class_numerosity, aes(x=Var1,y=Freq)) + geom_bar(stat="identity")

# too unbalanced!!

cut_dataset <- function(dataset){
  
  
  da <- dataset[dataset$Subtype=="LumA",]
  db <- dataset[dataset$Subtype=="LumB",]
  dh <- dataset[dataset$Subtype=="Her2",]
  ds <- dataset[dataset$Subtype=="Basal",]
  
  la<-length(da[,1])
  lb<-length(db[,1])
  lh<-length(dh[,1])
  ls<-length(ds[,1])
  
  mean <- (la+lb+lh+ls)/4
  
  if(la > mean) da <- da[sample(1:la,mean),]
  if(lb > mean) db <- db[sample(1:lb,mean),]
  if(lh > mean) dh <- dh[sample(1:lh,mean),]
  if(ls > mean) ds <- ds[sample(1:ls,mean),]
  
  rbind(da,db,dh,ds)
  
}

balanced_data <- cut_dataset(dataset)
  
class_numerosity = data.frame(table( balanced_data$Subtype))
ggplot(data=class_numerosity, aes(x=Var1,y=Freq)) + geom_bar(stat="identity")

# again?

balanced2 <- cut_dataset(balanced_data)
class_numerosity = data.frame(table( balanced2$Subtype))
ggplot(data=class_numerosity, aes(x=Var1,y=Freq)) + geom_bar(stat="identity")

# seem good

dataset <- balanced2
remove(balanced_data,balanced2,class_numerosity)


########################## HYPERPARAMETER TUNING

best_parameters <- NA

if(!file.exists("result/best_parameters.csv")){
  
  model_tuning <- boost_tree(mode="classification",
                             mtry=tune(),trees = tune(),
                             min_n=tune(),tree_depth=tune(),
                             learn_rate=tune(),loss_reduction=tune(),
                             sample_size=tune()) %>% set_engine("xgboost", objective = "multi:softprob", num_class=4)
  
  #let's multithread
  library(doParallel)
  all_cores = parallel::detectCores(logical = F)
  cl = makeCluster(5 ,setup_strategy="sequential")
  registerDoParallel(cl)
  
  
  tune_grid <- grid_regular(tree_depth(), trees(c(10,1000)),
                            min_n(),
                            loss_reduction(),
                            sample_size = sample_prop(),
                            finalize(mtry(), dataset),
                            learn_rate())
  
  tuning_result <- workflow() %>% add_formula(Subtype ~ .)  %>% add_model(model_tuning) %>% 
                            tune_grid(resamples = vfold_cv(dataset, 5),
                              grid = tune_grid, control=control_grid(parallel_over = "resamples"))
  
  
  #save bests
  best_parameters <- tuning_result %>% show_best("roc_auc")
  write.csv(best_parameters,"result/best_parameters.csv")
  
}else{
  
  best_parameters = read.csv("result/best_parameters.csv")
}



model_spec <-  boost_tree(mode="classification",
                          mtry= best_parameters$mtry[1], trees= best_parameters$trees[1],
                          min_n= best_parameters$min_n[1], tree_depth= best_parameters$tree_depth[1],
                          learn_rate= best_parameters$learn_rate[1], loss_reduction= best_parameters$loss_reduction[1],
                          sample_size= best_parameters$sample_size[1]) %>% set_engine("xgboost", objective = "multi:softprob", num_class=4)


################# VALIDATION

split <- initial_split(dataset,4/5)

model_eval <- model_spec %>% fit(Subtype ~ . , training(split))


prediction <- predict(model_eval,testing(split))
result <- cbind(predicted=prediction,actual=testing(split)$Subtype)

conf_mat <- conf_mat(result,truth=actual,estimate = .pred_class)
conf_mat <- conf_mat$table


# good!!

############### TRAINING

model <- model_spec %>% fit(Subtype ~ . , dataset)
saveRDS(model,"result/model_xgboost.rds")





