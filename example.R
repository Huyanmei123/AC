library(mclust)
source("./R/adaptive.R")


data<-read.csv("./goolam/goolam.csv",header = T,sep = ",", row.names = 1)
data<-as.matrix(data)
data<-log(data+1)


label<-read.csv("./goolam/goolam_label.csv", header = TRUE, row.names = 1)

#k: the number of the clusters for spectral clustering, if k is null, it will automatically get a value in the algorithm
#ncores: the number of workers 
#seed: spectral clustering
result <- adaptive(k = NULL, data = data, ncores = 4, seed = 0)

label <- as.vector(label$cell_type1)
final <- result$cluster

print(mclust::adjustedRandIndex(final,label))
