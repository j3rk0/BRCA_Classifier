library(gficf)

gficfNormalize <- function (matrix){
  gficf(matrix, cell_proportion_max = 1, cell_proportion_min = 0.01, storeRaw = TRUE, normalize = TRUE)
}

cluster_pca_tsne <- function (normalized){
  pca <- runPCA(normalized,50)
  tsne <- runReduction(pca,reduction = "tsne",seed = 0,nt=4)  #                                                     community detection
  clustcells(tsne, dist.method = "manhattan", nt = 4, k = 61 , community.algo = "louvian 2",seed=0,resolution = 0.5,n.start = 25,n.iter = 50)
}

cluster_pca_umap <- function (normalized){
  pca <- runPCA(normalized,50)
  umap <- runReduction(pca,reduction = "tumap",seed = 0,nt=4)
  clustcells(umap, dist.method = "manhattan", nt = 4, k = 61 , community.algo = "louvian 2",seed=0,resolution = 0.5,n.start = 25,n.iter = 50)
}

cluster_lsa_tsne <- function (normalized) {
  lsa <- runLSA(normalized,50)
  tsne <- runReduction(lsa,reduction = "tsne",seed = 0,nt=4)
  clustcells(tsne, dist.method = "manhattan", nt = 4, k = 61 , community.algo = "louvian 2",seed=0,resolution = 0.5,n.start = 25,n.iter = 50)
}

cluster_lsa_umap <- function (normalized){
  lsa <- runLSA(normalized,50)
  umap <- runReduction(lsa,reduction = "tumap",seed = 0,nt=4)
  clustcells(umap, dist.method = "manhattan", nt = 4, k = 61 , community.algo = "louvian 2",seed=0,resolution = 0.5,n.start = 25,n.iter = 50)
}

plot_cell_cluster <- function (clustered){
  plotCells(clustered,colorBy = "cluster")+ xlab("t-SNE1")+ylab("t-SNE2")
}

get_clust_matrix <- function (clust){
  barcodes <- colnames(clust$rawCounts)
  positions <- cells_indexes[match(barcodes,cells_indexes$barcode),]
  result <- cbind(positions,clust$embedded$cluster)
  colnames(result)<- c("barcode","x","y","cluster")
  result
}