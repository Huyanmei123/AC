library(doParallel)
library(foreach)



# Jaccard 相似度函数
jaccord <- function (x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  if(length(x) != length(y))
    stop("arguments must be vectors of the same length")
  tab <- table(x, y)
  if(all(dim(tab) == c(1, 1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  
  return(a / (a + b + c))
}

# 参数设置
a <- c(0.7, 0.75, 0.8, 0.85, 0.9)
b <- c(0.000000001, 0.000000005, 0.00000001, 0.00000005, 0.0000001, 0.0000005, 0.000001, 0.000005, 0.00001, 0.00005)
t <- c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)
dataSet <- c("campbell")

# 读取label
label <- read.csv("./data/clustering_result/campbell_label.csv", header = TRUE)
#######转为向量
label <- as.vector(label$cell_type)


# 设置并行
numCores <- 9
cl <- makeCluster(numCores)
registerDoParallel(cl)

# 仅对a并行
resultList <- foreach(tempa = a, .packages = c("Matrix", "aricode","cluster")) %dopar% {
  source("./R/Support.R")
  for (tempb in b) {
    for (tempt in t) {
      # 读取文件并处理
      result <- list()
      result$all <- lapply(1:9, function(index) {
        file <- paste("./data/clustering_result/campbell/campbell_fcom_", index, "_", tempa, "_", format(tempb, scientific = FALSE), "_", tempt, ".txt", sep = "")
        one_result_clustering_content <- read.table(file, header = FALSE, sep = "\t")
        one_result_clustering <- lapply(one_result_clustering_content, function(x) {
          as.integer(unlist(strsplit(x, " ")))
        })
        return(one_result_clustering)
      })
      print(paste0(""))
      
      # 设置seed保证重复性
      set.seed(123)
      final <- clustercom2(result)
      cluster <- final
      
      log_file <- file("clustering_test_log.txt", open = "a")
      writeLines(paste("length: ",), log_file)
      close(log_file)
      
    
      # 计算指标
      ari <- aricode::ARI(cluster, label)
      ji  <- jaccord(cluster, label)
      nmi <- aricode::NMI(cluster, label)
      
      # 日志记录
      log_file <- file("clustering_test_log.txt", open = "a")
      writeLines(paste("a:", tempa, "b:", tempb, "t:", tempt, "ari:", ari, "ji:", ji, "nmi:", nmi), log_file)
      close(log_file)
      
      # 将结果立即保存到文件
      write.table(data.frame(dataSet = "campbell", a = tempa, b = tempb, t = tempt, ari = ari, ji = ji, nmi = nmi),
                  file = "clustering_test_results.csv", 
                  sep = ",", 
                  col.names = !file.exists("clustering_test_results.csv"), 
                  row.names = FALSE, 
                  append = TRUE)
      gc()
    }
  }
}

# 停止并行集群
stopCluster(cl)
