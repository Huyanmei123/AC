library(mclust)
source("./R/adaptive.R")


data<-read.csv("./goolam/goolam.csv",header = T,sep = ",", row.names = 1)
data<-as.matrix(data)
data<-log(data+1)


label<-read.csv("./goolam/goolam_label.csv", header = TRUE, row.names = 1)
ncores <- 8
seed <- 0
result <- adaptive(k = NULL, data = data, ncores = ncores, seed = seed)

label <- as.vector(label$cell_type1)
final <- result$cluster

print(mclust::adjustedRandIndex(final,label))
