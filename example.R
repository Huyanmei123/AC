library(mclust)
source("./R/adaptive.R")

## load the goolam dataset
data<-read.csv("./goolam/goolam.csv",header = T,sep = ",", row.names = 1)
data<-as.matrix(data)

## Log transform the data
data<-log(data+1)

## load the labels of goolam label
label<-read.csv("./goolam/goolam_label.csv", header = TRUE, row.names = 1)

# k: the number of the clusters for spectral clustering, if k is null, it will automatically get a value in the algorithm
# ncores: the number of workers 
# seed: random seed
result <- adaptive(k = NULL, data = data, ncores = 4, seed = 0)

label <- as.vector(label$cell_type1)
final <- result$cluster

## evaluate
print(mclust::adjustedRandIndex(final,label))
