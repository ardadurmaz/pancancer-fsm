#!/usr/bin/env Rscript

ConnComponent <- function(ids = NULL, graph = NULL){
  if(length(ids) < 2){
    return(NULL)
  }
  idx <- colnames(graph) %in% ids
  local.graph <- graph[idx, idx]
  visited <- setNames(numeric(ncol(local.graph)), colnames(local.graph))
  
  compNum <- 0
  ConnComponentRes <- list()
  for(i in 1:ncol(local.graph)){
    if(visited[i] == 0){
      visited[i] <- 1
      compNum <- compNum + 1
      queue <- c(names(visited)[i])
      ConnComponent <- c()
      while(length(queue) > 0){
        vertex <- queue[length(queue)]
        queue <- queue[-length(queue)]
        ConnComponent <- c(ConnComponent,vertex)
        
        idx <- which(local.graph[,which(colnames(local.graph) == vertex)] == 1)
        if(length(idx) > 0){
          for(k in 1:length(idx)){
            if(visited[idx[k]] == 0){
              visited[idx[k]] <- 1
              queue <- c(names(visited)[idx[k]], queue)
            }
          }
        }
      }
      ConnComponentRes[[compNum]] <- ConnComponent
    }
  }
  ConnComponentRes <- ConnComponentRes[unlist(lapply(ConnComponentRes, length)) >= 8]
  return(ConnComponentRes)
}


suppressMessages(library(Matrix))
suppressMessages(library(parallel))
suppressMessages(library(optparse))


option_list <- list(make_option(c("-r", "--profile_results"), action = "store", type = "character", help = "Network Profile Results", default = NULL),
                    make_option(c("-n", "--network_data"), action = "store", type = "character", help = "Network Data", default = NULL),
                    make_option(c("-t", "--num_threads"), action = "store", type = "integer", help = "Number of Threads", default = 2),
                    make_option(c("-o", "--output_prefix"), action = "store", type = "character", help = "Prefix of Output Folder", default = "RWRGraphs"),
                    make_option(c("-s", "--support"), action = "store", type = "integer", help = "Minimum Support", default = 8))

opt <- parse_args(OptionParser(option_list = option_list))
output.dir <- sprintf("%s/%s", getwd(), opt$output_prefix)


## Pre-Process ##
if(is.null(opt$profile_results) || is.null(opt$network_data)){
  warning('!! Provide random walk results and network data')
  quit(save = "no")
}else{
  tryCatch({
    dir.create(output.dir)
  }, warning = function(w){
    cat(sprintf("!! Warning [%s]", w))
    quit(status = 1, save = "no")
  },error = function(err){
    cat(sprintf("!! Error [%s]", e))
    quit(status = 1, save = "no")
  })
}


message("...reading data")
rwr.res <- readRDS(opt$profile_results)
graph <- readRDS(opt$network_data)

sig.genes <- apply(rwr.res, 2, function(x){
  names(x)[x < 0.1]
})
rm(rwr.res)
gc()

idx <- sapply(sig.genes, function(x){length(x) == 0})
if(sum(idx) > 0){
  warning(sprintf("!! %d samples with no significant features", sum(idx)))
  sig.genes <- sig.genes[!idx]
}

message("...extracting connected components")
cl <- makeCluster(opt$num_threads)
clusterExport(cl, varlist = c('ConnComponent', 'graph'))
cnn <- parSapply(cl, sig.genes, simplify = FALSE, function(ids){
  ConnComponent(ids, graph)
})
stopCluster(cl)
gc()

message("...sorting network")
idx <- order(colnames(graph))
graph <- graph[idx, idx]
edge.idx <- which(Matrix::triu(graph), arr.ind = TRUE, useNames = FALSE)

message("...extracting edges")
cl <- makeCluster(opt$num_threads)
clusterExport(cl, varlist = c('cnn', 'graph', 'edge.idx'))
sig.edges <- parSapply(cl, cnn, simplify = FALSE, function(cnn.l){
  local.ids <- sapply(cnn.l, simplify = FALSE, function(genes.l){
    which(colnames(graph)[edge.idx[,1]] %in% genes.l & colnames(graph)[edge.idx[,2]] %in% genes.l)
  })
  return(Reduce(union, local.ids))
})
stopCluster(cl)
sig.edges <- sig.edges[!sapply(sig.edges, is.null)]

message("...calculating edge support levels")
edge.counts <- numeric(nrow(edge.idx))
for(i in 1:length(sig.edges)){
  local.sig.edges <- sig.edges[[i]]
  if(is.null(local.sig.edges))
    next
  edge.counts[local.sig.edges] <- edge.counts[local.sig.edges] + 1
}
edge.counts.nSupport <- which(edge.counts >= opt$support)

message("..writing sample specific graphs")
cl <- makeCluster(opt$num_threads)
clusterExport(cl, varlist = c('sig.edges', 'graph', 'edge.idx', 'edge.counts.nSupport', 'output.dir'))
written.res <- parSapply(cl, 1:length(sig.edges), function(s){
  sample.edges <- intersect(edge.counts.nSupport, sig.edges[[s]])
  if(length(sample.edges) > 3){
    sample.graph <- sapply(1:length(sig.edges), simplify = FALSE, function(i){
      local.edges <- intersect(sig.edges[[i]], sample.edges)
      if(length(local.edges) == 0 || is.null(local.edges)){
        return(NULL)
      }
      local.genes <- union(unique(edge.idx[local.edges, 1]), 
                           unique(edge.idx[local.edges, 2]))
      names(local.genes) <- colnames(graph)[local.genes]
      off.set <- (i-1)*ncol(graph)
      
      vertices <- sapply(1:length(local.genes), simplify = FALSE, function(g.i){
        val <- sprintf("vertex\t{\"id\":%d,\"label\":\"%s\"}\n", 
                       local.genes[g.i] + off.set, 
                       names(local.genes)[g.i])
      })
      vertices <- do.call('c', lapply(vertices, function(x){return(x)}))
      
      edges <- sapply(local.edges, simplify = FALSE, function(e.i){
        temp <- edge.idx[e.i,]
        temp.genes <- colnames(graph)[temp]
        val.1 <- sprintf("edge\t{\"src\":%d,\"targ\":%d,\"label\":\"%s\"}\n", 
                         local.genes[which(names(local.genes) == temp.genes[1])] + off.set, 
                         local.genes[which(names(local.genes) == temp.genes[2])] + off.set, 
                         paste(temp.genes, collapse = "."))
        
        return(val.1)
      })
      edges <- do.call('c', lapply(edges, function(x){return(x)}))
      return(c(vertices,edges))
    })
    sample.graph <- do.call('c', lapply(sample.graph, function(x){return(x)}))
    cat(sample.graph, file = sprintf('%s/%s_graph.tsv', output.dir, names(sig.edges)[s]), sep = '', append = FALSE)
  }
  return('Done')
})
stopCluster(cl)

quit(save = "no", status = 0)
