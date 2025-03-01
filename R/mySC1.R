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
 
  norm_X <- sqrt(rowSums(X^2))
  norm_Y <- sqrt(rowSums(Y^2))
  
  
  dot_product <- tcrossprod(X, Y)
  
 
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

    h <- num[i, ]
    
 
    hj <- seq(1, nrow(num))
    

    hj <- hj[hj != i]
    
   
    h <- h[hj]
    
 
    sorted_indices <- order(h, decreasing = TRUE)
    h <- h[sorted_indices]
    hj <- hj[sorted_indices]
    
   
    min_index <- 0
    if (is.null(kmax)) {
      m <- mean(h)
      bool_tensor <- h < m
   
      endindex <- which.max(bool_tensor)
      h <- h[1:endindex]
      hj <- hj[1:endindex]
      
      # 
      hh <- diff(h)
      hhh <- diff(hh)
      
    
      min_index <- which.min(hhh)
    } else {
      min_index <- kmax
    }
    
 
    selected_hj <- hj[1:min_index]
    selected_h <- h[1:min_index]
    
   
    D <- c(D, list(selected_h))
    C <- c(C, list(selected_hj))
  }
  DC = list(D, C)
 
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
  m = length(DC[[1]])  
  edges <- matrix(ncol = 3, nrow = 0)  
  for(i in 1:m) {
    for(j in 1:length(DC[[2]][[i]])) {
      edge1 <- c(i, DC[[2]][[i]][[j]], DC[[1]][[i]][[j]])
      edges <- rbind(edges, edge1)  
    }
  }
  
  
  g <- NULL
  if(ncol(edges) == 3) {
    edges_df <- as.data.frame(edges)
    colnames(edges_df) <- c("from", "to", "weight")
    g <- graph_from_data_frame(as.data.frame(edges_df), directed = FALSE)
  } else {
    stop("Edge list format error: edges should have 3 columns (source, target, weight).")
  }
  
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
  #data <- as.matrix(data) 
  da = apply(data,1, paste, collapse="#")
  indUnique = which(!duplicated(da))
  indAll = match(da, da[indUnique])
  
  data2 = data
  data  = data[indUnique, ]
  n <- nrow(data)
  
 
  
  
  if(is.null(gmax)){
    if(!is.null(centers)) gmax = centers - 1L
    else gmax = 1L
  }
  test=TRUE
 
  
  DC.tmp = mydist(data)

  DC = list(DC.tmp[[1]], DC.tmp[[2]])
  W <- DC[[1]]
 
  wi <- c()
  for(i in 1:length(W))
  {
    wi <- c(wi,W[[i]][length(W[[i]])])
  }

  SC <- lapply(DC[[2]], function(sublist) lapply(sublist, function(x) rep(1, length(x))))
  for(i in 1:length(W))
  {
    for(j in 1:length(W[[i]]))
    {
      SC[[i]][[j]] <- wi[DC[[2]][[i]][[j]]][[1]] * wi[[i]]
    }
  }

  for(i in 1:length(SC))
  {
    for(j in 1:length(SC[[i]]))
    {
      W[[i]][[j]] <- (W[[i]][[j]]*W[[i]][[j]]) / SC[[i]][[j]]
    }
  }

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
