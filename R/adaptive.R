
source('./R/Support.R')
library(Matrix)
library(foreach)
library(infotheo)
library(torch)
library(Matrix)
library(data.table)
library(scDHA)
library(doParallel)
library(foreach)



adaptive <- function(k = NULL, data, method = "adaptive", sparse = FALSE,n = 5e3, ncores = 10L, gen_fil = TRUE, sample.prob = NULL, seed = NULL) 
{
  gc()
  print("programming start")
  res <- adaptive.basic(k = k, method = method, data=data, n = n, ncores = ncores, gen_fil = gen_fil, sample.prob = sample.prob, seed = seed)
  res
}

adaptive.basic <- function(k = NULL, data, method = "adaptive", n = 5e3, ncores = 10L, gen_fil = TRUE, sample.prob = NULL, seed = NULL)
{
  
    latent <- scDHA(data = data, do.clus = FALSE, seed=1)
    latent <- latent$all.latent
   
 
    result <- list()
    cl <- parallel::makeCluster(ncores, outfile = "/dev/null")
    doParallel::registerDoParallel(cl, cores = ncores)
    parallel::clusterEvalQ(cl,{
      library(scDHA)
    })
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
    g.en <- latent[[which.max(sapply(result$all, function(x) adjustedRandIndex(x,final)))]]
    final <- as.numeric(factor(final))
    
    
    
    list( cluster = final,
          latent = g.en,
          all.latent = latent,
          all.res = result$all)
    # keep.genes = keep
  
  
}
