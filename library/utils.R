normalize_seq_rawcount <- function(matrix){
  
  library(edgeR)
  
  matrix <- t(matrix)
  matrix<- cpm(calcNormFactors(DGEList(counts=matrix),normalized.lib.sizes = T))
  matrix<-log10(matrix+1)
  t(matrix)
}


calculate_roc <- function(class_prob){
  
  res <- data.frame(x=NA,y=NA)
  res <- res[-1,]
  
  for (i in 0:10) {
    
    treshold <- (1/10) * i
    
    tp <- 0
    fp <- 0
    tn <- 0
    fn <- 0
    
    for (j in 1:nrow(class_prob)) {
      if(class_prob$prob[j] > treshold){
        if(class_prob$actual[j] == 1 ) tp <- tp +1
        else fp <- fp + 1
      }else{
        if(class_prob$actual[j] == 0 ) tn <- tn +1
        else fn <- fn + 1 
      }
    }
    
    spec <- tn / (tn+fp)
    sens <- tp / (tp+fn)
    
    res <- rbind(res,c(1-spec,sens))
  }
  
  colnames(res) <- c("x","y")
  res
}


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


evalutate_model_performance <- function(model_spec,dataset,folds,prefix,n_folds){
  
  #final rate matrix
  
  rm <- matrix(rep(0,20),nrow = 4,ncol = 4)
  rownames(rm)<-c("LumA","LumB","Basal","Her2")#,"Normal")
  colnames(rm)<-c("tpr","tnr","fpr","fnr")
  
  #for each fold
  
  for (i in 1:n_folds) {
    
    #take indexes
    training_indexes <- unlist(folds[i])
    
    #train model
    model <- model_spec %>% fit(Subtype ~ ., data= dataset[ training_indexes,])
    #predict
    predictions <- predict(model,dataset[-training_indexes,])
    
    #create a dataset predicted/actual class
    truth <- dataset$Subtype[-training_indexes]
    data<- as.data.frame(cbind(predictions,truth))
    
    #build confusion matrix
    conf_mat <- conf_mat(data=data,truth=truth,estimate = .pred_class)$table
    
    #export confusion matrix
    path <- paste0(prefix,"/f",i,"/confmtx.csv")
    dir.create(path)
    write.table(conf_mat,file=path)
    
    #fill rate matrix with counts
    rm["LumA","tpr"]<-rm["LumA","tpr"]+ conf_mat["LumA","LumA"]
    rm["LumA","tnr"]<-rm["LumA","tnr"]+ sum(conf_mat) - sum(conf_mat["LumA",])-sum( conf_mat[,"LumA"]) + conf_mat["LumA","LumA"]
    rm["LumA","fpr"]<-rm["LumA","fpr"]+ sum(conf_mat["LumA",]) - conf_mat["LumA","LumA"]
    rm["LumA","fnr"]<-rm["LumA","fnr"]+ sum(conf_mat[,"LumA"]) - conf_mat["LumA","LumA"]
    rm["LumB","tpr"]<-rm["LumB","tpr"]+ conf_mat["LumB","LumB"]
    rm["LumB","tnr"]<-rm["LumB","tnr"]+ sum(conf_mat) - sum(conf_mat[,"LumB"])-sum( conf_mat["LumB",])+ conf_mat["LumB","LumB"]
    rm["LumB","fpr"]<-rm["LumB","fpr"]+ sum(conf_mat["LumB",]) - conf_mat["LumB","LumB"]  
    rm["LumB","fnr"]<-rm["LumB","fnr"]+ sum(conf_mat[,"LumB"])-conf_mat["LumB","LumB"]
    rm["Basal","tpr"]<-rm["Basal","tpr"]+ conf_mat["Basal","Basal"]
    rm["Basal","tnr"]<-rm["Basal","tnr"]+sum(conf_mat)-sum(conf_mat["Basal",])-sum(conf_mat[,"Basal"])+ conf_mat["Basal","Basal"]
    rm["Basal","fpr"]<-rm["Basal","fpr"]+sum(conf_mat["Basal",])-conf_mat["Basal","Basal"]
    rm["Basal","fnr"]<-rm["Basal","fnr"]+ sum(conf_mat[,"Basal"])-conf_mat["Basal","Basal"]
    rm["Her2","tpr"]<-rm["Her2","tpr"]+ conf_mat["Her2","Her2"]
    rm["Her2","tnr"]<-rm["Her2","tnr"]+sum(conf_mat)-sum(conf_mat["Her2",])-sum(conf_mat[,"Her2"])+ conf_mat["Her2","Her2"]
    rm["Her2","fpr"]<-rm["Her2","fpr"]+sum(conf_mat["Her2",])-conf_mat["Her2","Her2"]
    rm["Her2","fnr"]<-rm["Her2","fnr"]+ sum(conf_mat[,"Her2"])-conf_mat["Her2","Her2"]
  }
  
  #switch from counts to rating
  
  rm["LumA","tpr"] <- rm["LumA","tpr"] / (rm["LumA","tpr"] + rm["LumA","fnr"])
  rm["LumA","fpr"] <- rm["LumA","fpr"] / (rm["LumA","fpr"] + rm["LumA","tnr"])
  rm["LumA","tnr"] <- rm["LumA","tnr"] / (rm["LumA","tnr"] + rm["LumA","fpr"])
  rm["LumA","fnr"] <- rm["LumA","fnr"] / (rm["LumA","fnr"] + rm["LumA","tpr"])
  
  rm["LumB","tpr"] <- rm["LumB","tpr"] / (rm["LumB","tpr"] + rm["LumB","fnr"])
  rm["LumB","fpr"] <- rm["LumB","fpr"] / (rm["LumB","fpr"] + rm["LumB","tnr"])
  rm["LumB","tnr"] <- rm["LumB","tnr"] / (rm["LumB","tnr"] + rm["LumB","fpr"])
  rm["LumB","fnr"] <- rm["LumB","fnr"] / (rm["LumB","fnr"] + rm["LumB","tpr"])
  
  rm["Basal","tpr"] <- rm["Basal","tpr"] / (rm["Basal","tpr"] + rm["Basal","fnr"])
  rm["Basal","fpr"] <- rm["Basal","fpr"] / (rm["Basal","fpr"] + rm["Basal","tnr"])
  rm["Basal","tnr"] <- rm["Basal","tnr"] / (rm["Basal","tnr"] + rm["Basal","fpr"])
  rm["Basal","fnr"] <- rm["Basal","fnr"] / (rm["Basal","fnr"] + rm["Basal","tpr"])
  
  rm["Her2","tpr"] <- rm["Her2","tpr"] / (rm["Her2","tpr"] + rm["Her2","fnr"])
  rm["Her2","fpr"] <- rm["Her2","fpr"] / (rm["Her2","fpr"] + rm["Her2","tnr"])
  rm["Her2","tnr"] <- rm["Her2","tnr"] / (rm["Her2","tnr"] + rm["Her2","fpr"])
  rm["Her2","fnr"] <- rm["Her2","fnr"] / (rm["Her2","fnr"] + rm["Her2","tpr"])
  
  #export rating matrix
  
  path<-paste0(prefix,"/result.csv")
  
  
  write.csv(as.data.frame(rm),path)
  
}


plot_dataset_distribution <- function(dataset){
  ggplot(data=dataset,aes(Subtype)) + geom_bar()
}


plot_compare_roc <-  function(r1,r2,par){
  
  roc_1<- as.data.frame( cbind( prob= r1[,par] , actual= as.numeric(r1$actual_class==par) ))
  roc_1 <- calculate_roc(roc_1)
  
  roc_2 <- as.data.frame( cbind( prob= r2[,par], actual= as.numeric(r2$actual_class==par) ))
  roc_2 <- calculate_roc( roc_2 )
  
  toPlot <- roc_1 %>% mutate(params= "curve1") %>% bind_rows(roc_2 %>% mutate(params="curve2"))
  ggplot(toPlot,aes(x=x,y=y,color=params)) + geom_line() + xlab("1-specificity") + ylab("sensitivity") + ggtitle(paste(par," roc curve"))
}
