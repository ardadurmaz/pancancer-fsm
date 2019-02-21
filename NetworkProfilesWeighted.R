#!/usr/bin/env Rscript

suppressMessages(library(Matrix))
suppressMessages(library(parallel))
suppressMessages(library(optparse))

option_list <- list(make_option(c("-m", "--mutation_data"), action = "store", type = "character", help = "Mutation Matrix", default = NULL),
                    make_option(c("-e", '--expression_data'), action = "store", type = "character", help = "Expression Matrix", default = NULL),
                    make_option(c("-n", "--network_data"), action = "store", type = "character", help = "Network Data", default = NULL),
                    make_option(c("-t", "--num_threads"), action = "store", type = "integer", help = "Number of Threads", default = 2),
                    make_option(c("-o", "--output_prefix"), action = "store", type = "character", help = "Output Prefix", default = "NetworkProfile"))

opt <- parse_args(OptionParser(option_list = option_list))

if(is.null(opt$mutation_data) || 
   is.null(opt$network_data) ||
   is.null(opt$expression_data)){
  warning("!Provide data files")
  quit(save = "no")
}

message("...loading data")
mut.mat <- readRDS(opt$mutation_data)
expr.mat <- readRDS(opt$expression_data)
ppi <- readRDS(opt$network_data)

message("...processing data")
if(length(intersect(rownames(mut.mat),
                    intersect(rownames(expr.mat),
                              rownames(ppi)))) == 0){
  warning("!Feature names do not match")
  quit(save = "no")
}

## Filter Samples ##
samples.common <- intersect(colnames(mut.mat), colnames(expr.mat))
mut.mat <- mut.mat[,colnames(mut.mat) %in% samples.common]
expr.mat <- expr.mat[,colnames(expr.mat) %in% samples.common]

## Filter Features ##
idx <- colnames(ppi) %in% na.omit(intersect(rownames(mut.mat), rownames(expr.mat)))
ppi <- ppi[idx, idx]
idx <- colSums(ppi) == 0
ppi <- ppi[!idx,!idx]

mut.mat <- mut.mat[match(rownames(ppi),table=rownames(mut.mat)),]
mut.mat <- mut.mat[,colSums(mut.mat) > 8]
expr.mat <- expr.mat[match(rownames(ppi),table=rownames(expr.mat)),]
expr.mat <- expr.mat[,match(colnames(mut.mat),table=colnames(expr.mat))]
expr.mat <- t(scale(t(expr.mat)))
gc()

network.prof <- sapply(1:ncol(mut.mat), function(i){
  cat(sprintf('\rProcessing Sample: %d',i))
  
  ## Weighted Network ##
  edge.idx <- which(ppi, arr.ind = TRUE, useNames = FALSE)
  edge.idx <- cbind(edge.idx, apply(edge.idx, 1, function(x){
    sqrt(expr.mat[x[1],i]^2 + expr.mat[x[2],i]^2)
  }))
  edge.idx[is.na(edge.idx[,3]),3] <- 0.1 ## Small dysregulation
  ppi.w <- sparseMatrix(i = edge.idx[,1],
                        j = edge.idx[,2],
                        x = edge.idx[,3],
                        dims = c(nrow(ppi), nrow(ppi)),
                        dimnames = list(rownames(ppi), rownames(ppi)),
                        use.last.ij = TRUE)
  ppi.w <- ppi.w %*% Matrix::Diagonal(x = Matrix::colSums(ppi.w)^-1)
  
  ## Original ##
  p0 <- Matrix(mut.mat[,i], ncol = 1, nrow = nrow(ppi.w))
  p0 <- p0 / sum(p0)
  pt <- Matrix(1/nrow(ppi.w), ncol = 1, nrow = nrow(ppi.w))
  
  delta <- 1
  count <- 1
  r <- 0.75
  while(delta > 1e-16 && count < 100){
    px <- (1-r) * ppi.w %*% pt + r * p0
    delta <- sum(abs(px - pt))
    count <- count + 1
    pt <- px
  }
  res.orig <- pt
  
  ## Random ##
  nSeed <- sum(mut.mat[,i])
  cl <- makeCluster(opt$num_threads)
  clusterExport(cl, varlist = c('ppi.w', 'r', 'nSeed'), envir = environment())
  res.random <- parSapply(cl, 1:1000, function(k){
    require(Matrix)
    p0 <- Matrix(0, ncol = 1, nrow = nrow(ppi.w))
    p0[sample(1:nrow(ppi.w), size = nSeed, replace = FALSE), 1] <- 1/nSeed
    pt <- Matrix(1/nrow(ppi.w), ncol = 1, nrow = nrow(ppi.w))
    delta <- 1
    count <- 1
    while(delta > 1e-16 && count < 100){
      px <- (1-r) * ppi.w %*% pt + r * p0
      delta <- sum(abs(px - pt))
      count <- count + 1
      pt <- px
    }
    return(as.vector(pt))
  })
  stopCluster(cl)
  
  ## Significance ##
  res.random <- log10(res.random)
  mean.vals <- apply(res.random, 1, function(x){mean(x, na.rm = TRUE)})
  sd.vals <- apply(res.random, 1, function(x){sd(x, na.rm = TRUE)})
  p.vals <- pnorm(log10(as.vector(res.orig)),
                  mean = mean.vals,
                  sd = sd.vals, 
                  lower.tail = FALSE)
  p.vals <- p.adjust(p.vals, method = 'fdr')
  return(setNames(p.vals, nm = rownames(ppi.w)))
})
colnames(network.prof) <- colnames(mut.mat)

out.file <- sprintf("%s/%s_NetworkProfile.rds", getwd(), opt$output_prefix)
message("...saving results to file ", out.file)
saveRDS(network.prof, file = out.file)
quit(save = "no")