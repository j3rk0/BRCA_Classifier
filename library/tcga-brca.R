library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(GenomicDataCommons)


load_brca_experiment_data <- function (){
  query <- GDCquery("TCGA-BRCA",   #tcga breast carcinoma
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts",
                    legacy = FALSE
  )
  
  GDCdownload(query)
  GDCprepare(query)
}

#return rawCOunt from experiment
get_brca_genes_count <- function (Experiment){
  assay <- assay(Experiment)
  genes <- rowRanges(Experiment)
  genes <- ranges(genes)
  genes <- names(genes)
  row.names(assay) <- genes
  as.data.frame(t(assay))
}

#return clinicaldata from experiment
get_brca_clinic_data <- function (experiment){
  data <- colData(experiment)
  
}

get_brca_dataframe <- function (experiment) {
  
  gc <- get_brca_genes_count(experiment)
  cd <- get_brca_clinic_data(experiment)
  res <- as.data.frame(cbind(cd,gc))
  
}

#return a dataframe with pam50 genes count and cancer subtype
get_brca_pam50_dataset <- function(experiment){
  
  source("library/geneOperation.R")
  source("library/pam50.R")
  
  cd <- get_brca_clinic_data(experiment)
  gc <- filter_gene_count_matrix(get_brca_genes_count(experiment),get_pam50_ids())
  
  subtype <- cd$paper_BRCA_Subtype_PAM50
  as.data.frame(cbind(subtype,gc))
  
}
