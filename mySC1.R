#' @importFrom igraph arpack decompose graph
#' @importFrom stats cor dnorm qnorm

fast.table <- function (data)
{
  if(!is.data.frame(data))
    data = as.data.frame(data, stringsAsFactors = FALSE)
  da = do.call("paste", c(data, sep = "\r"))
  ind = !duplicated(da)
  levels = da[ind]
  cat <- factor(da,levels = levels)
  nl <- length(levels(cat))
  bin <- (as.integer(cat) - 1)
  pd <- nl
  bin <- bin[!is.na(bin)]
  if (length(bin)) bin <- bin + 1
  y <- tabulate(bin, pd)
  result=list(index = bin, weights = y, data = data[ind,])
  result
}

my_cos <- function(X, Y) {
  # 将 X 和 Y 转换为 Torch 张量
  norm_X <- sqrt(rowSums(X^2))
  norm_Y <- sqrt(rowSums(Y^2))
  
  # 计算X和Y的内积矩阵
  dot_product <- tcrossprod(X, Y)
  
  # 使用outer来归一化内积结果，得到余弦相似度矩阵
  m <- dot_product / outer(norm_X, norm_Y)
  
  return(m)
}

createFolds <- function(vec, k)
{
  rand.vec <- sample.int(length(vec), replace = FALSE, n = length(vec))
  tmp <- round(seq(1, max(vec), length.out = k+1))
  res <- list()
  for (i in 1:k) {
    res[[i]] <- rand.vec[tmp[i]:tmp[i+1]]
  }
  res
}

origin_mydist <- function(data, k = 20, distance = 2)
{

  m <- dim(data)[1]
  q <- dim(data)[2]
  D <- matrix(nrow = m, ncol = k)
  C <- matrix(nrow = m, ncol = k)
  folds <- createFolds(1:m, k = ceiling(m/1000))
  
  for (idx in folds) {
    tmp <- data[idx,]
    dis.tmp <- 1 - cor(t(tmp), t(data))
    
    for (i in 1:length(idx)) {
      tmp1 <- order(dis.tmp[i,])
      tmp1 <- tmp1[2:(k+1)]
      
      D[idx[i],] <- dis.tmp[i,tmp1]
      C[idx[i],] <- tmp1
    }
  }
  
  list(D,C)
  
}


mydist <- function(data)
{
  kmax = NULL
  num = data
  num = my_cos(num, num)
  D <- list()
  C <- list()
  
  for (i in 1:nrow(num)) {
    # if (i %% 1000 == 0) {
    #   cat("Processing:", i, "\n")
    # }
    
    # 获取第i行数据
    h <- num[i, ]
    
    # 创建一个从1到nrow(num)的序列，表示所有的行索引
    hj <- seq(1, nrow(num))
    
    # 排除当前行
    hj <- hj[hj != i]
    
    # 根据排除后的索引来选择数据
    h <- h[hj]
    
    # 按照降序排序
    sorted_indices <- order(h, decreasing = TRUE)
    h <- h[sorted_indices]
    hj <- hj[sorted_indices]
    
    # 计算均值并进行过滤
    min_index <- 0
    if (is.null(kmax)) {
      m <- mean(h)
      bool_tensor <- h < m
      #endindex <- 250  # 获取第一个小于均值的最大索引
      endindex <- which.max(bool_tensor)
      h <- h[1:endindex]
      hj <- hj[1:endindex]
      
      # 计算一阶和二阶导数
      hh <- diff(h)
      hhh <- diff(hh)
      
      # 获取二阶导数的最小值索引
      min_index <- which.min(hhh)
    } else {
      min_index <- kmax
    }
    
    # 选择前min_index个元素
    selected_hj <- hj[1:min_index]
    selected_h <- h[1:min_index]
    
    # 将选择的结果添加到D和C
    D <- c(D, list(selected_h))
    C <- c(C, list(selected_hj))
  }
  DC = list(D, C)
  cat("mydist\n")
  return(DC)
}


getClosest = function(X, Y)
{
  m = nrow(Y)
  n = nrow(X)
  res = matrix(0, n, m)
  for(i in 1:m){
    tmp = (X-rep(Y[i,], each=n))**2
    res[,i] = rowSums(tmp)
  }
  apply(res, 2, which.min)
}


Laplacian <- function(DC, k, normalize="none") {
  m = length(DC[[1]])  # 获取DC[[1]]的长度
  edges <- matrix(ncol = 3, nrow = 0)  # 创建一个空矩阵，3列用于保存边和权重
  print(paste0("Laplacian graph ", "construct"))
  for(i in 1:m) {
    for(j in 1:length(DC[[2]][[i]])) {
      edge1 <- c(i, DC[[2]][[i]][[j]], DC[[1]][[i]][[j]])
      edges <- rbind(edges, edge1)  # 将边信息按行追加到矩阵中
    }
  }
  
  # 确保edges是正确的格式
  g <- NULL
  if(ncol(edges) == 3) {
    edges_df <- as.data.frame(edges)
    colnames(edges_df) <- c("from", "to", "weight")
    g <- graph_from_data_frame(as.data.frame(edges_df), directed = FALSE)
  } else {
    stop("Edge list format error: edges should have 3 columns (source, target, weight).")
  }
  
  # 计算拉普拉斯矩阵
  #edges_list <- list()  # 存储边的列表
  # print(paste0("Laplacian graph ", "construct"))
  # # 添加边
  # for (i in 1:m) {
  #   for (j in 1:length(DC[[2]][[i]])) {
  #     edge1 <- c(i, DC[[2]][[i]][[j]], DC[[1]][[i]][[j]])
  #     edge2 <- c(DC[[2]][[i]][[j]], i, DC[[1]][[i]][[j]])
  #     
  #     # 将每条边转化为唯一的字符标识符
  #     edge1_key <- paste(edge1, collapse = "_")
  #     edge2_key <- paste(edge2, collapse = "_")
  #     
  #     # 检查是否已存在，若不存在，则添加
  #     if (!exists(edge1_key, envir = edge_set)) {
  #       assign(edge1_key, TRUE, envir = edge_set)  # 标记已存在
  #       edges_list <- append(edges_list, list(edge1))  # 添加到边列表
  #     }
  #     
  #     if (!exists(edge2_key, envir = edge_set)) {
  #       assign(edge2_key, TRUE, envir = edge_set)
  #       edges_list <- append(edges_list, list(edge2))
  #     }
  #   }
  # }
  # matrix <- do.call(rbind, lapply(edges_list, function(x) unlist(x)))
  # print(paste0("Laplacian graph ", "finish"))
  print(paste0("Laplacian graph ", "finish"))
  ##########获得g的source边
  source = get.edgelist(g)
  weights = E(g)$weight
  i = c(as.numeric(source[,1]), as.numeric(source[,2]))
  j = c(as.numeric(source[,2]), as.numeric(source[,1]))
  x = c(weights, weights)  # 对称图，权重重复添加一次
  matrix = cbind(i,j,x)
  result <- sparseMatrix(i = matrix[,1], j = matrix[,2], x = matrix[,3], dims = c(m, m))
  
  # Step 4: Compute the degree vector D
  D <- rowSums(result)
  
  # Step 5: Apply normalization
  if (normalize == "none") {
    return(Diagonal(x = D) - result)
  }
  
  if (normalize == "symmetric") {
    TMP <- Diagonal(x = 1 / sqrt(D))
    result <- TMP %*% result %*% TMP
    return(Diagonal(m) - result)
  }
  
  if (normalize == "random-walk") {
    return(Diagonal(m) - Diagonal(x = 1 / D) %*% result)
  }
  
  return(result)
}


AUC = function(y)
{
  l = length(y)
  x = 0:(l-1)
  y = y - y[1]
  res = numeric(0)
  
  for(i in 1:l){
    A = 0
    A = y[i]*(i-1)/2
    B = y[i] * (l-i)
    C = (y[l] - y[i]) *  (l-i) / 2
    res[i] = A+B+C
  }
  res
}




specClust <- function (data, centers=NULL, nn = 7, method = "symmetric", gmax=NULL, ...)
{
 
  call = match.call()
  if(is.data.frame(data)) data = as.matrix(data)
  # unique data points
  da = apply(data,1, paste, collapse="#")
  indUnique = which(!duplicated(da))
  indAll = match(da, da[indUnique])
  
  data2 = data
  data  = data[indUnique, ]
  n <- nrow(data)
  
  #data = scale(data, FALSE, TRUE)
  
  
  if(is.null(gmax)){
    if(!is.null(centers)) gmax = centers - 1L
    else gmax = 1L
  }
  test=TRUE
  # DC.tmp = origin_mydist(data, 30)
  # while(test){
  #   if(nn > ncol(DC.tmp[[1]])) DC.tmp = origin_mydist(data, nn*2)
  #   
  #   DC = list(DC.tmp[[1]][,1:nn], DC.tmp[[2]][,1:nn])
  #   sif <- rbind(1:n, as.vector(DC[[2]]))
  #   g <- graph(sif, directed=FALSE)
  #   g <- decompose(g, min.vertices=4)
  #   if (length(g) > 1) {
  #     #warning("graph not connected")
  #     if(length(g)>=gmax) nn = nn+2
  #     else test=FALSE
  #   }
  #   else test=FALSE
  # }
  
  DC.tmp = mydist(data)
  # cat(paste0("nn=",nn,"\n"),file = "./outputdata/logger.txt",append = T)
  DC = list(DC.tmp[[1]], DC.tmp[[2]])
  W <- DC[[1]]
  # #######选取W每一个list的最后一个元素
  wi <- c()
  for(i in 1:length(W))
  {
    wi <- c(wi,W[[i]][length(W[[i]])])
  }
  ##(W)    DC[[2]]
  SC <- lapply(DC[[2]], function(sublist) lapply(sublist, function(x) rep(1, length(x))))
  for(i in 1:length(W))
  {
    for(j in 1:length(W[[i]]))
    {
      SC[[i]][[j]] <- wi[DC[[2]][[i]][[j]]][[1]] * wi[[i]]
    }
  }
  # SC[] <-  wi[DC[[2]]] * wi
  for(i in 1:length(SC))
  {
    for(j in 1:length(SC[[i]]))
    {
      W[[i]][[j]] <- (W[[i]][[j]]*W[[i]][[j]]) / SC[[i]][[j]]
    }
  }
  # W = W^2 / SC
  for(i in 1:length(W))
  {
    alpha = 1/(2*(length(W[[i]])+1))
    qua = abs(qnorm(alpha))
    for(j in 1:length(W[[i]]))
    {
      W[[i]][[j]] = W[[i]][[j]]*qua
      W[[i]][[j]] = dnorm(W[[i]][[j]], sd = 1)
    }
  }
  
  DC[[1]] = W
  
  L = Laplacian(DC, nn, method)
 
  
  f <- function(x, extra) as.vector(extra %*% x)
  
  if(is.null(centers))kmax = 25
  else kmax = max(centers)
  
  U <- arpack(f, extra = L, options = list(n = n, which = "SM",
                                           nev = kmax, ncv = 5 * kmax, mode=1), sym = TRUE)
 
  ind <- order(U[[1]])
  U[[2]] = U[[2]][indAll, ind]
  U[[1]] = U[[1]][ind]
  if (is.null(centers)) {
    tmp = which.max(diff(U[[1]]))+1
    centers = which.min(AUC(U[[1]][1:tmp]))
  }
  if(method == "symmetric"){
    rs = sqrt(rowSums(U[[2]]^2))
    U[[2]] =  U[[2]]/rs
  }
  result = kmeans(U[[2]], centers = centers, nstart = 50, iter.max = 5000, ...)

  archeType = getClosest(U[[2]][indAll, ], result$centers)
  result$eigenvalue = U[[1]]
  result$eigenvector = U[[2]]
  result$data = data2
  result$indAll = indAll
  result$indUnique = indUnique
  result$L = L
  result$archetype = archeType
  result$call = call
  class(result) = c("specClust", "kmeans")
  result
}