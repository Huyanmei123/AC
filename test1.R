
source('./R/Support.R')
library(Matrix)
library(foreach)
library(infotheo)
library(aricode)
library(torch)
library(Matrix)
library(data.table)


k = NULL
seed = 2
method = "scDHA"
do.clus = TRUE
ncores <- 9L


jaccord <- function (x, y)
{
  x <- as.vector(x)
  y <- as.vector(y)
  if(length(x) != length(y))
    stop("arguments must be vectors of the same length")
  tab <- table(x,y)
  if(all(dim(tab)==c(1,1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c

  return(a/(a+b+c))
}

adaptive <- function(k = NULL, method = "scDHA", sparse = FALSE, dataset,n = 5e3, ncores = 10L, gen_fil = TRUE, do.clus = TRUE, sample.prob = NULL, seed = NULL) 
{
  gc()
  res <- adaptive.basic(k = k, method = method, n = n, ncores = ncores, gen_fil = gen_fil, do.clus = do.clus, sample.prob = sample.prob, seed = seed,dataset=dataset)
  res
}

adaptive.basic <- function(k = NULL, method = "scDHA", n = 5e3, ncores = 10L, gen_fil = TRUE, do.clus = TRUE, sample.prob = NULL, seed = NULL,dataset)
{
 
  load("./outputdata/xin.RData")
  latent <- result$all.latent
  
  # 循环读取文件并合并
  # for (index in 1:9) {
  #   #temp <- read.csv(paste0("./outputdata/myoutput/", dataset, "/", dataset, "_data_", index, ".csv"), header = TRUE, sep = ",")
  #   
  #   # 将读取的数据追加到 'latent' 中
  #   #latent <-  c(latent, list(as.matrix(temp)))
  #  
  # }
  # 
  if (do.clus) {
    result <- list()
    cl <- parallel::makeCluster(ncores, outfile = "/dev/null")
    doParallel::registerDoParallel(cl, cores = ncores)
    parallel::clusterEvalQ(cl,{
      library(scDHA)
    })
    print("第三个多进程")
    result$all <- lapply(latent, function(x) {
      source('./R/Support.R')
      library(Matrix)
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
      
      if (method == "louvain") {
        set.seed(seed)
        cluster <- as.numeric(factor(clus.louvain(x)))
      } 
      else {
        set.seed(seed)
        if (nrow(x) >= 50e3) {
          cluster <- clus.big(x, k = k, n = 5e3, nmax = 15)
          cluster <- as.numeric(factor(cluster))
        } else {
          if (nrow(x) > 5e3) {
            cluster <- clus.big(x, k = k)
          } else {
            cluster <- clus(x, k = k)
          }
        }
      }
      
      return(cluster)
    })
    set.seed(seed)
    final <- clustercom2(result)
    print(final)
    g.en <- latent[[which.max(sapply(result$all, function(x) adjustedRandIndex(x,final)))]]
    final <- as.numeric(factor(final))
    
    list( cluster = final,
          latent = g.en,
          all.latent = latent,
          all.res = result$all)
    # keep.genes = keep
  } else {
    list( all.latent = latent,
          keep.genes = keep)
  }
  
}

dataset = "xin"
sink("logfile.txt")
result <- adaptive(k = NULL, method = method, n = 5e3, ncores = ncores, gen_fil = TRUE, do.clus = do.clus, sample.prob = NULL, seed = seed,dataset = dataset)


# 读取label
label <- read.csv(paste0("./data/clustering_result/",dataset,"_label.csv"), header = TRUE)
#######转为向量
label <- as.vector(label$cell_type1)
final <- result$cluster
ari <- aricode::ARI(final,label)
ji  <- jaccord(final,label)
nmi <- aricode::NMI(final,label)
print(ari)
print(ji)
print(nmi)






# dataSet<-c("yan","goolam","deng","pollen","patel","wang","darmanis","camp-brain",
#            "usoskin","kolodziejczyk","camp-liver","xin","baron-mouse","muraro","segerstolpe")
dataSet<-c("yan")
# 
testResult <- data.frame(dataSet,ari=0,ji=0,nmi=0)
for(data in dataSet)
{
  if(is.null(colnames(data))) keep.genes <- seq(ncol(data)) else keep.genes <- colnames(data)
  # 加载数据集
  #a <- load(paste0("./data/",data,".RData"))
  #######打开.csv文件，并且取SVM_clusterID这一列
  label <- read.csv("./data/clustering_result/chen_label.csv", header = TRUE, row.names = 1)
  label <- label$SVM_clusterID
  
  print(label)
  cat(paste0("data=",data,"\n"),file = "./outputdata/logger.txt",append = T)
  
  # 读取9个降维矩阵
  
  
 
  cluster <- result$cluster

  ari <- aricode::ARI(cluster,label)
  ji  <- jaccord(cluster,label)
  nmi <- aricode::NMI(cluster,label)
  print(ari)
  testResult[testResult$dataSet==data,"ari"] <- ari
  testResult[testResult$dataSet==data,"ji"] <- ji
  testResult[testResult$dataSet==data,"nmi"] <- nmi
  gc()
  
}

write.csv(testResult,"./outputdata/testResult.csv",row.names = FALSE)



latent <- result$all.latent
rm(result)
gc()
source("./R/Analysis.R")
library(viridis)
result = NULL
result$all <- lapply(1:9, function(index) {
  cluster <- read.table(paste0("./data/clustering_result/","patel","/","patel","_",index,".txt"),header = F,sep="")
  #########一行转数组
  cluster <- as.vector(t(cluster))
  
  cluster
})
final <- clustercom2(result)
latent1 <- latent[[which.max(sapply(result$all, function(x) adjustedRandIndex(x,final)))]]
result$pred <- uwot::umap(latent1, n_threads = ncores)
library(RColorBrewer)
plot(result$pred, 
     col =  rainbow(length(unique(result$pred)))[factor(result$pred)],  # Viridis color palette
     xlab = "KNN-1", 
     ylab = "KNN-2", 
     xlim = range(-0.9, 0.9), 
     ylim = range(-0.9, 0.9), 
     pch = 19, 
     cex = 1.2,  # Point size
     xaxt = "n",  # Disable default x-axis
     yaxt = "n"   # Disable default y-axis
)
axis(1, at = seq(-0.9, 0.9, length.out = 5), cex.axis = 1)
axis(2, at = seq(-0.9, 0.9, length.out = 5), cex.axis = 1)
#' data <- log2(data + 1)
#' if(torch::torch_is_installed()) #scDHA need libtorch installed
#' {
#'   #Generate clustering result, the input matrix has rows as samples and columns as genes
#'   result <- scDHA(data, ncores = 2, seed = 1)
#'   #Generate 2D representation, the input is the output from scDHA function
#'   result <- scDHA.vis(result, ncores = 2, seed = 1)
#'   #Plot the representation of the dataset, different colors represent different cell types
#'   plot(result$pred, col=factor(label), xlab = "scDHA1", ylab = "scDHA2")
#' }
#' }
#'
label <- read.csv(paste0("./data/clustering_result/","patel","_label.csv"), header = TRUE)
label <- as.vector(label$cell_type1)
library(scDHA)
data <- patel@assays$data$logcounts
data <- as.matrix(data)
# resultCluster<-scDHA(data=data,seed = 1)
resultCluster<- adaptive.basic(k = NULL, method = "scDHA", n = 5e3, ncores = 10L, gen_fil = TRUE, do.clus = TRUE, sample.prob = NULL, seed = 1, dataset = "patel")
result <- scDHA.vis(resultCluster, ncores = 10L, seed = 1)
plot(result$pred, col=viridis(length(unique(label)))[factor(label)], xlab = "AC-1", ylab = "AC-2",xlim=range(-20,20),ylim=range(-20,20),pch = 19, cex = 1.2, xaxt = "n", yaxt = "n")
axis(1, at = seq(min(range(-20,20)), max(range(-20,20)), length.out = 5), cex.axis = 1.2)
axis(2, at = seq(min(range(-20,20)), max(range(-20,20)), length.out = 5), cex.axis = 1.2)


library(Seurat)
library(scDHA)

for(index in 1:9){
  write.csv(result$all.latent[[index]],paste0("./outputdata/myoutput/patel/","patel_data_",index,".csv"),row.names = F)
}



result = NULL
result$all <- lapply(1:9, function(index) {
  cluster <- read.table(paste0("./data/clustering_result/","patel","/","patel","_",index,"_knn12.txt"),header = F,sep="")
  #########一行转数组
  cluster <- as.vector(t(cluster))
  
  cluster
})
final <- clustercom2(result)
write.table(final,"patel_knn12.txt",sep = "\t",row.names = F,col.names = F)
