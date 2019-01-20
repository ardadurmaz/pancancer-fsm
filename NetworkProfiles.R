#!/usr/bin/env Rscript

prob.vec <- function(graph = NULL, r = 0.75, nThread = 2){
  require(Matrix)
  require(parallel)
  
  ppi.n <- graph %*% Matrix::Diagonal(x = Matrix::colSums(graph)^-1)
  cl <- makeCluster(nThread)
  clusterExport(cl, varlist = c('ppi.n', 'r'), envir = environment())
  rwr.res <- parSapply(cl, 1:ncol(ppi.n), function(i){
    require(Matrix)
    p0 <- Matrix(0, ncol = 1, nrow = nrow(ppi.n))
    p0[i, 1] <- 1
    pT <- Matrix(1/nrow(ppi.n), ncol = 1, nrow = nrow(ppi.n))
    delta <- 1
    count <- 1
    while(delta > 1e-16 && count < 100){
      pX <- (1-r) * ppi.n %*% pT + r * p0
      delta <- sum(abs(pX - pT))
      count <- count + 1
      pT <- pX
    }
    return(as.vector(pT))
  })
  stopCluster(cl)
  colnames(rwr.res) <- colnames(ppi.n)
  rownames(rwr.res) <- colnames(ppi.n)
  return(rwr.res)
}

rwr <- function(rwr.prof = NULL, seed = NULL, nThread = 2){
  require(Matrix)
  require(parallel)
  
  cl <- makeCluster(nThread)
  clusterExport(cl, varlist = c('rwr.prof'), envir = environment())
  rwr.res <- parApply(cl, seed, 2, function(x){
    local.seed <- names(x)[x != 0]
    pt.orig <- Matrix::rowMeans(rwr.prof[,which(colnames(rwr.prof) %in% local.seed)])
    pt.random <- sapply(1:1000, function(i){
      return(as.vector(Matrix::rowMeans(rwr.prof[,sample(1:ncol(rwr.prof), size = length(local.seed), replace = FALSE)])))
    })
    pt.random <- log(pt.random)
    mean.stats <- apply(pt.random, 1, mean)
    sd.stats <- apply(pt.random, 1, sd)
    test.res <- setNames(log(as.vector(pt.orig)), nm = rownames(pt.orig))
    p.vals <- pnorm(test.res, mean = mean.stats, sd = sd.stats, lower.tail = FALSE)
    p.vals.adj <- p.adjust(p.vals, method = 'bonferroni')
    return(p.vals.adj)
  })
  stopCluster(cl)
  gc()
  rownames(rwr.res) <- rownames(seed)
  return(rwr.res)
}



suppressMessages(library(Matrix))
suppressMessages(library(parallel))
suppressMessages(library(optparse))

option_list <- list(make_option(c("-m", "--mutation_data"), action = "store", type = "character", help = "Mutation Matrix", default = NULL),
                    make_option(c("-n", "--network_data"), action = "store", type = "character", help = "Network Data", default = NULL),
                    make_option(c("-t", "--num_threads"), action = "store", type = "integer", help = "Number of Threads", default = 2),
                    make_option(c("-o", "--output_prefix"), action = "store", type = "character", help = "Output Prefix", default = "NetworkProfile"))

opt <- parse_args(OptionParser(option_list = option_list))

if(is.null(opt$mutation_data) || is.null(opt$network_data)){
  warning("!Provide mutation and network data")
  quit(save = "no")
}

message("...loading mutation and network data")
mut.mat <- readRDS(opt$mutation_data)
ppi <- readRDS(opt$network_data)

message("...processing data")
if(sum(rownames(mut.mat) %in% rownames(ppi)) == 0){
  warning("!Feature names do not match")
  quit(save = "no")
}

mut.mat.ft <- apply(mut.mat, 2, function(x){
  temp <- setNames(numeric(nrow(ppi)), nm = rownames(ppi))
  temp[which(names(temp) %in% names(x)[x == 1])] <- 1
  return(temp)
})
idx <- colSums(mut.mat.ft) == 0
if(sum(idx) > 0){
  warning(sprintf("!Removing %d samples which do not have features represented", sum(idx)))
  mut.mat.ft <- mut.mat.ft[,!idx]
}

message("...calculating probability vectors")
rwr.prof <- prob.vec(graph = ppi, r = 0.6, nThread = opt$num_threads)

message("...generating network profiles")
network.prof <- rwr(rwr.prof = rwr.prof, seed = mut.mat.ft, nThread = opt$num_threads)

out.file <- sprintf("%s/%s_NetworkProfile.rds", getwd(), opt$output_prefix)
message("...saving results to file ", out.file)
saveRDS(network.prof, file = out.file)
quit(save = "no")