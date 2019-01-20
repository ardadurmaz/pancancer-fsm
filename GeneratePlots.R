library(Matrix)
library(stringi)
library(parallel)

setwd('~/pancancer_fsm/')
cl <- makeCluster(20)
fsm.res <- parSapply(cl,
                     list.files(path = '/home/arda/pancancer_fsm/results/combined_s20',
                                pattern = 'names.txt', 
                                full.names = TRUE, 
                                recursive = TRUE,
                                include.dirs = TRUE),
                     simplify = FALSE,
                     function(fName){
                       require(stringi)
                       if(file.size(fName) > 0){
                         tryCatch({
                           local.fsg.data <- read.table(fName, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
                         },
                         warning = function(c){
                           warning("Issue reading file")
                           return(NULL)
                         },
                         error = function(c){
                           return(NULL)
                         })
                         sample.name <- as.character(stri_match(fName, regex = 'TCGA-.{2}-.{4}-\\d{2}', mode = 'first'))
                         return(data.frame('sample.id' = rep(sample.name, 
                                                             nrow(local.fsg.data)),
                                           'fsg' = local.fsg.data,
                                           stringsAsFactors = FALSE))
                       }else{
                         return(NULL)
                       }
                     })
stopCluster(cl)

fsm.res <- do.call('rbind',
                   lapply(fsm.res,
                          function(x){
                            return(x)
                          }))
colnames(fsm.res) <- c('sample.id', 'fsg')
temp <- split(fsm.res$fsg, fsm.res$sample.id)

cl <- makeCluster(20)
clusterExport(cl, varlist = c('temp'))
fsm.genes.pp <- parSapply(cl, 1:length(temp), simplify = FALSE, function(x){
  require(stringi)
  local.s <- stri_match_all(temp[[x]], regex = "\\d\\:([^\\(]+?)\\]")
  local.g <- lapply(local.s, function(k){
    local.local.g <- paste(sort(unique(unlist(strsplit(k[,2],split = '.',fixed = TRUE)))), collapse = ',')
  })
  local.g <- unique(do.call('c', lapply(local.g, function(k){return(k)})))
  return(data.frame('sample.id' = rep(names(temp)[x], length(local.g)),
                    'fsg' = local.g,
                    stringsAsFactors = FALSE))
})
stopCluster(cl)
fsm.genes.pp.ft <- do.call('rbind', lapply(fsm.genes.pp, function(x){return(x)}))

sig.fsg <- fsm.genes.pp.ft$fsg[which(as.numeric(table(fsm.genes.pp.ft$fsg)) > 8)]
sig.fsg.pp <- sapply(sig.fsg, simplify = FALSE, function(local.fsg){
  fsm.genes.pp.ft$sample.id[which(fsm.genes.pp.ft$fsg == local.fsg)]
})
all.sample.ids <- unique(fsm.genes.pp.ft$sample.id)
sig.fsg.pp.mat <- matrix(0, ncol = length(all.sample.ids), nrow = length(sig.fsg))

for(i in 1:length(sig.fsg.pp)){
  sig.fsg.pp.mat[i, which(all.sample.ids %in% sig.fsg.pp[[i]])] <- 1
}

idx <- colSums(sig.fsg.pp.mat) == 0
sig.fsg.pp.mat <- sig.fsg.pp.mat[,!idx]


require(fastcluster)
require(parallelDist)
require(pheatmap)

hc.dist <- parDist(t(sig.fsg.pp.mat), method = 'manhattan')
hr.dist <- parDist(sig.fsg.pp.mat, method = 'manhattan')

hc <- hclust(hc.dist, method = 'average')
hr <- hclust(hr.dist, method = 'average')

pheatmap(sig.fsg.pp.mat[hr$order,hc$order],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 12,
         show_rownames = FALSE,
         show_colnames = FALSE,
         filename = '~/FsmPlots/PatientFSGCluster.pdf',
         width = 12,
         height = 10)

fsg.genes.unique <- unique(do.call("c", lapply(stri_match_all(fsm.res$fsg,
                                                              regex = "\\d\\:([^\\(]+?)\\]"),
                                               function(x){
                                                 local.s <- x[,2]
                                                 return(paste(sort(unique(unlist(strsplit(local.s,
                                                                                          split='.',
                                                                                          fixed=TRUE)))),
                                                              collapse=','))
                                                 })))
## Graph Cluster ##
fsm.edge <- unique(do.call("c", lapply(stri_match_all(fsm.res$fsg, 
                                                      regex = "\\d\\:([^\\(]+?)\\]"), 
                                       function(x){x[,2]})))
fsm.edge.ft <- do.call('rbind', strsplit(fsm.edge, split = '.', fixed = TRUE))
fsm.edge.ft <- as.data.frame(fsm.edge.ft, stringsAsFactors = FALSE)

require(Matrix)
require(igraph)

colnames(fsm.edge.ft) <- c('gene.a', 'gene.b')
all.ids <- union(fsm.edge.ft$gene.a,
                 fsm.edge.ft$gene.b)
fsm.edge.ft <- rbind(fsm.edge.ft,
                     data.frame('gene.a' = fsm.edge.ft$gene.b,
                                'gene.b' = fsm.edge.ft$gene.a,
                                stringsAsFactors = FALSE))
row.idx <- match(fsm.edge.ft$gene.a, table = all.ids)
col.idx <- match(fsm.edge.ft$gene.b, table = all.ids)
adj <- sparseMatrix(i = row.idx,
                    j = col.idx,
                    x = 1,
                    dims = c(length(all.ids), length(all.ids)),
                    dimnames = list(all.ids, all.ids),
                    use.last.ij = TRUE)
graph <- graph_from_adjacency_matrix(adj, mode = 'undirected')
graph.cluster <- cluster_louvain(graph)
cluster.data <- data.frame('symbol' = names(membership(graph.cluster)),
                           'cluster.id' = as.vector(membership(graph.cluster)),
                           stringsAsFactors = FALSE)
saveRDS(cluster.data, file = '~/pancancer_fsm/results/ClusterDataSupport20.rds')

library(ggplot2)
cluster.size <- data.frame('size' = as.vector(sizes(graph.cluster)),
                           'cluster' = paste('Cluster-', 1:length(as.vector(sizes(graph.cluster))), sep = ''))
p <- ggplot(cluster.size,
       aes(x = reorder(cluster.size$cluster, X = cluster.size$size), y = size)) +
  geom_bar(stat = 'identity',
           fill = 'grey50') + 
  theme_minimal() +
  xlab("Cluster") +
  ylab("Size") +
  ggtitle("Cluster Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12))
ggsave(p,
       file = '~/FsmPlots/ClusterSizeSupport20.pdf',
       width = 10,
       height = 8)

## Cluster Enrichment ##
reactome.data <- read.table('~/pancancer_fsm/data/Ensembl2Reactome_All_Levels.txt',
                            header = FALSE,
                            sep = '\t',
                            stringsAsFactors = FALSE,
                            quote = '',
                            comment.char = '')
reactome.data <- subset(reactome.data, reactome.data$V6 == 'Homo sapiens')

library(biomaRt)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(mart = ensembl,
                 attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id",
                 values = unique(reactome.data$V1))
reactome.data$V1 <- mapping$hgnc_symbol[match(reactome.data$V1, table = mapping$ensembl_gene_id)]
reactome.list <- split(reactome.data$V1, reactome.data$V2)
reactome.list <- lapply(reactome.list, function(x){unique(na.omit(x))})

ppi <- readRDS("~/pancancer_fsm/data/stringdb.hsa.hgnc.medconf.adj.rds")
reactome.list <- lapply(reactome.list, function(x){
  x[x %in% colnames(ppi)]
})
reactome.list <- reactome.list[sapply(reactome.list, length) >= 8]

graph.cluster.list <- split(cluster.data$symbol, 
                            cluster.data$cluster.id)
names(graph.cluster.list) <- paste('Cluster-', names(graph.cluster.list), sep = '')
graph.cluster.list <- graph.cluster.list[sapply(graph.cluster.list, length) >= 8]
graph.cluster.list <- graph.cluster.list[sapply(graph.cluster.list, length) < 400]

univ <- colnames(ppi)
cl <- makeCluster(20)
clusterExport(cl, varlist = c('univ', 'reactome.list'))
enrichRes <- parSapply(cl, 
                       graph.cluster.list, 
                       simplify = FALSE,
                       function(x){
                         local.enr <- sapply(reactome.list,
                                             simplify = TRUE,
                                             function(t){
                                               test.mat <- matrix(c(sum(x %in% t),
                                                                    sum(!(x %in% t)),
                                                                    sum(univ %in% t),
                                                                    sum(!(univ %in% t))),
                                                                  byrow = TRUE, ncol = 2)
                                               return(fisher.test(test.mat, alternative = 'greater')$p.value)
                                             })
                         return(p.adjust(local.enr, method = 'bonferroni'))
                       })
stopCluster(cl)
enrichRes.ft <- do.call('cbind', lapply(enrichRes, function(x){return(x)}))

p.thr <- 0.01
enrichRes.ft <- enrichRes.ft[rowSums(enrichRes.ft < p.thr) > 0,]
enrichRes.ft <- -1 * log10(enrichRes.ft)
rownames(enrichRes.ft) <- reactome.data$V4[match(rownames(enrichRes.ft), table = reactome.data$V2)]

## Reactome Hierarchy ##
reactome.h <- read.table('data/ReactomePathwaysRelation.txt',
                         header = FALSE,
                         sep = '\t',
                         stringsAsFactors = FALSE)
reactome.h <- reactome.h[reactome.h$V1 %in% names(reactome.list) & reactome.h$V2 %in% names(reactome.list),]
root.nodes <- setdiff(reactome.h$V1, reactome.h$V2)
require(igraph)
graph <- graph_from_edgelist(as.matrix(reactome.h), directed = TRUE)

short.path <- sapply(setdiff(V(graph)$name, root.nodes), simplify = TRUE, function(x){
  local.path <- shortest_paths(graph, from = x, to = root.nodes, mode = 'in')
  local.path <- local.path[!is.null(local.path)]
  local.path.ft <- local.path$vpath
  local.dist <- setNames(sapply(local.path.ft, function(k){length(k$name)}), nm = root.nodes)
  local.dist <- local.dist[local.dist > 0]
  return(names(local.dist)[which.min(local.dist)])
})

reactome.mapping <- data.frame('Node' = names(short.path),
                               'Root.Node' = short.path,
                               stringsAsFactors = FALSE)
reactome.mapping <- rbind(reactome.mapping,
                          data.frame('Node' = reactome.mapping$Root.Node,
                                     'Root.Node' = reactome.mapping$Root.Node,
                                     stringsAsFactors = FALSE))

require(reshape2)
temp.data <- melt(enrichRes.ft)
colnames(temp.data) <- c('Pathway', 'Cluster', 'Enrichment (Log10 Scale)')
temp.data <- subset(temp.data, temp.data$`Enrichment (Log10 Scale)` > 2)
write.table(temp.data, file = '~/pancancer_fsm/results/ClusterPathwayEnrichmentTableSupport20.tsv', col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE, sep = '\t')


# ## Plot ##
top.pathways <- apply(enrichRes.ft, 2, function(x){order(x, decreasing = TRUE)[1:3]})
top.pathways <- unique(as.vector(top.pathways))


library(pheatmap)
library(RColorBrewer)
breaksList = seq(0, 200, by = 1)
colors <- colorRampPalette(brewer.pal(n = 9, "Reds"))(length(breaksList)-1)
colors <- colors[50:200]
breaksList <- seq(0, 200, by = 200 / length(colors))
enrichRes.ft.filt <- enrichRes.ft[top.pathways,]


temp.mapp <- reactome.data$V2[match(rownames(enrichRes.ft.filt), table = reactome.data$V4)]
temp.temp.mapp <- myReactomeHier$root[match(temp.mapp, table = myReactomeHier$child)]
temp.temp.temp.mapp <- reactome.data$V4[match(temp.temp.mapp, table = reactome.data$V2)]

col_annotation <- data.frame("Root.Pathway" = temp.temp.temp.mapp,
                             row.names = rownames(enrichRes.ft.filt))
colors.manual <- c('aliceblue','antiquewhite3','antiquewhite4',
                   'aquamarine', 'aquamarine3', 'aquamarine4',
                   'azure', 'azure4', 'beige', 
                   'blueviolet', 'burlywood', 'burlywood4',
                   'cadetblue', 'cadetblue4', 'chartreuse', 
                   'chartreuse4', 'chocolate', 'cornflowerblue', 
                   'cyan', 'cyan4', 'darkgoldenrod1', 
                   'darkorchid4', 'darkred','hotpink1', 
                   'darksalmon', 'darkslategrey', 'greenyellow',
                   'lightpink4', 'seashell4', 'blue4',
                   'yellow', 'tan', 'slateblue4')


require(fastcluster)
hc <- hclust(dist(t(enrichRes.ft.filt)), method = 'average')
hr <- hclust(dist(enrichRes.ft.filt), method = 'average')


row.ordering <- order(col_annotation$Root.Pathway)
col.ordering <- sapply(row.ordering, simplify = FALSE, function(i){
  which(enrichRes.ft.filt[i, ] > 2)
})
col.ordering <- unique(do.call("c", lapply(col.ordering, function(x){return(x)})))



enrichRes.ft.filt <- enrichRes.ft.filt[row.ordering,
                                       col.ordering]

enrichRes.ft.filt <- t(enrichRes.ft.filt)

pheatmap(enrichRes.ft.filt, 
         fontsize = 12,
         annotation_col = col_annotation,
         annotation_colors = list('Root.Pathway' = setNames(colors.manual[1:length(unique(col_annotation$Root.Pathway))],
                                                            unique(col_annotation$Root.Pathway))),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_col = 4,
         fontsize_row = 6,  
         scale = "none",
         border_color = 'gray95',
         display_numbers = FALSE,
         show_colnames = FALSE,
         breaks = breaksList, 
         color = c("white",colors), 
         main = 'Reactome Pathway Enrichment',
         filename = '~/FsmPlots/PathwayEnrichmentSupport20v3.pdf',
         width = 12,
         height = 8)



## HotNet2 Enrichment ##
hotnet2 <- read.table('~/pancancer_fsm/data/HotNet2.tsv',
                      header = TRUE,
                      sep = '\t',
                      stringsAsFactors = FALSE,
                      fill = TRUE)
hotnet2 <- sapply(1:ncol(hotnet2), simplify = FALSE, function(x){
  local.res <- trimws(as.character(hotnet2[,x][hotnet2[,x] != '']))
  return(local.res[local.res %in% univ])
})
names(hotnet2) <- paste('Subnetwork', 1:length(hotnet2), sep = '')


## Cluster HotNet2 Enrichment ##
cl <- makeCluster(20)
clusterExport(cl, varlist = c('hotnet2', 'univ'))
enrichResClusterHotnet2 <- parSapply(cl, graph.cluster.list, simplify = FALSE, function(x){
  local.enr <- sapply(hotnet2, simplify = TRUE, function(t){
    test.mat <- matrix(c(sum(x %in% t),
                         sum(!x %in% t),
                         sum(univ %in% t),
                         sum(!univ %in% t)),
                       byrow = TRUE, ncol = 2)
    return(fisher.test(test.mat, alternative = 'greater')$p.value)
  })
  return(p.adjust(local.enr, method = 'bonferroni'))
})
stopCluster(cl)
enrichResClusterHotnet2.ft <- do.call('rbind', lapply(enrichResClusterHotnet2, function(x){return(x)}))
enrichResClusterHotnet2.ft <- enrichResClusterHotnet2.ft[,colSums(enrichResClusterHotnet2.ft < p.thr) > 0]
enrichResClusterHotnet2.ft <- enrichResClusterHotnet2.ft[rowSums(enrichResClusterHotnet2.ft < p.thr) > 0,]
enrichResClusterHotnet2.ft <- melt(enrichResClusterHotnet2.ft)
colnames(enrichResClusterHotnet2.ft) <- c('Cluster', 'Hotnet2 Subnetwork', 'p-value')
enrichResClusterHotnet2.ft <- subset(enrichResClusterHotnet2.ft, enrichResClusterHotnet2.ft$`p-value` < p.thr)
write.table(enrichResClusterHotnet2.ft, 
            file = '~/pancancer_fsm/results/ClusterHotnet2Enrichment.tsv',
            col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, quote = FALSE)

## Patient HotNet2 Enrichment ##
fsm.pp.gene <- sapply(split(fsm.res$fsg, fsm.res$sample.id), 
                      simplify = FALSE, 
                      function(x){
  local.genes <- unique(do.call("c",
                                lapply(stri_match_all(x, 
                                                      regex = "\\d\\:([^\\(]+?)\\]"), 
                                       function(m){
                                         paste(sort(unlist(strsplit(m[,2],split = '.',fixed = TRUE))), collapse = ',')
                                       })))
  return(local.genes)
})

cl <- makeCluster(20)
clusterExport(cl, varlist = c('hotnet2', 'univ'))
enrichRes.pp <- parSapply(cl, fsm.pp.gene, simplify = FALSE, function(x){
  local.enr <- sapply(x, simplify = TRUE, USE.NAMES = FALSE, function(x.local){
    x.local <- unlist(strsplit(x.local, split = ',', fixed = TRUE))
    local.local.enr <- sapply(hotnet2, simplify = TRUE, function(t){
        test.mat <- matrix(c(sum(x.local %in% t),
                             sum(!x.local %in% t),
                             sum(univ %in% t),
                             sum(!univ %in% t)),
                           byrow = TRUE, ncol = 2)
        return(fisher.test(test.mat, alternative = 'greater')$p.value)
    })
    return(p.adjust(local.local.enr, method = 'bonferroni'))
  })
  return(apply(local.enr, 1, function(x){min(x, na.rm = TRUE)}))
})
stopCluster(cl)

enrichRes.pp.ft <- do.call('rbind',
                           lapply(enrichRes.pp, function(x){
                             return(x)
                           }))
enrichRes.pp.ft <- as.matrix(t(enrichRes.pp.ft))
enrichRes.pp.ft <- -1 * log10(enrichRes.pp.ft)
temp <- melt(enrichRes.pp.ft)
colnames(temp) <- c('HotNet2 Subnetwork', 'Sample ID', 'Enrichment (Log Scale)')
temp <- subset(temp, temp$`Enrichment (Log Scale)` > 2)
write.table(temp, file = '~/pancancer_fsm/results/PatientHotNet2Enrichment.tsv', col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE, sep = '\t')

## Patient Cluster Enrichment ##
library(reshape2)

cl <- makeCluster(20)
clusterExport(cl, varlist = c('graph.cluster.list', 'univ'))
enrichRes.pp <- parSapply(cl, fsm.pp.gene, simplify = FALSE, function(x){
  local.enr <- sapply(x, simplify = TRUE, USE.NAMES = FALSE, function(x.local){
    x.local <- unlist(strsplit(x.local, split = ',', fixed = TRUE))
    local.local.enr <- sapply(graph.cluster.list, simplify = TRUE, function(t){
      test.mat <- matrix(c(sum(x.local %in% t),
                           sum(!(x.local %in% t)),
                           sum(univ %in% t),
                           sum(!(univ %in% t))),
                         byrow = TRUE, ncol = 2)
      return(fisher.test(test.mat, alternative = 'greater')$p.value)  
    })
    return(p.adjust(local.local.enr, method = 'bonferroni'))
  })
  return(apply(local.enr, 1, function(x){min(x, na.rm = TRUE)}))
})
stopCluster(cl)

enrichRes.pp.ft <- do.call('rbind',
                           lapply(enrichRes.pp, function(x){
                             return(x)
                           }))
enrichRes.pp.ft <- as.matrix(t(enrichRes.pp.ft))
enrichRes.pp.ft <- 1 - enrichRes.pp.ft
enrichRes.pp.ft[enrichRes.pp.ft > 0.99] <- 1
enrichRes.pp.ft[enrichRes.pp.ft <= 0.99] <- 0



pheno.data <- read.table('~/pancancer_fsm/data/TCGA_phenotype_denseDataOnlyDownload.tsv',
                         header = TRUE,
                         sep = '\t',
                         comment.char = '',
                         stringsAsFactors = FALSE)
col.annotation <- data.frame('Disease' = pheno.data$X_primary_disease[match(colnames(enrichRes.pp.ft), 
                                                                            table = pheno.data$sample)],
                             row.names = colnames(enrichRes.pp.ft))


colors.manual <- c('aliceblue','antiquewhite3','antiquewhite4',
                   'aquamarine', 'aquamarine3', 'aquamarine4',
                   'azure', 'azure4', 'beige', 
                   'blueviolet', 'burlywood', 'burlywood4',
                   'cadetblue', 'cadetblue4', 'chartreuse', 
                   'chartreuse4', 'chocolate', 'cornflowerblue', 
                   'cyan', 'cyan4', 'darkgoldenrod1', 
                   'darkorchid4', 'darkred','hotpink1', 
                   'darksalmon', 'darkslategrey', 'greenyellow',
                   'lightpink4', 'seashell4', 'blue4',
                   'yellow', 'tan', 'slateblue4')

breaksList = seq(0, 1, by = 0.5)
colors <- c('white', 'tomato')

require(fastcluster)
hc <- hclust(dist(t(enrichRes.pp.ft), method = 'manhattan'), method = 'average')
hr <- hclust(dist(enrichRes.pp.ft, method = 'manhattan'), method = 'average')

pheatmap(enrichRes.pp.ft[hr$order, hc$order],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = col.annotation,
         scale = "none",
         fontsize = 12,
         fontsize_col = 2,
         fontsize_row = 8,
         legend = FALSE,
         border_color = "gray95",
         breaks = breaksList,
         color = colors, 
         show_colnames = FALSE, 
         main = 'Patient Cluster Enrichment',
         annotation_colors = list("Disease" = setNames(colors.manual[1:length(unique(col.annotation$Disease))], 
                                                       unique(col.annotation$Disease))),
         filename = '~/FsmPlots/PatientClusterEnrichmentSupport20.pdf',
         width = 12,
         height = 10)





enrichRes.pp.ft <- melt(enrichRes.pp.ft)
enrichRes.pp.ft$phenotype <- pheno.data$X_primary_disease[match(enrichRes.pp.ft$Var2,
                                                                table = pheno.data$sample)]
colnames(enrichRes.pp.ft) <- c('cluster.id', 'sample.id', 'enr', 'disease')

enrichRes.pp.ft <- subset(enrichRes.pp.ft, enrichRes.pp.ft$enr == 1)
rownames(enrichRes.pp.ft) <- 1:nrow(enrichRes.pp.ft)

sample <- sapply(unique(enrichRes.pp.ft$cluster.id), simplify = FALSE, function(c.id){
  local.data <- subset(enrichRes.pp.ft, enrichRes.pp.ft$cluster.id == c.id)
  counts <- table(local.data$disease)
  return(data.frame('cluster.id' = paste('Cluster-', rep(c.id, length(counts)), sep = ''),
                    'disease' = names(counts),
                    'count' = as.vector(counts),
                    stringsAsFactors = FALSE))
})
sample.ft <- do.call('rbind', lapply(sample, function(x){return(x)}))


sample.ft$cluster.id <- factor(sample.ft$cluster.id, levels = paste('Cluster-', sort(unique(as.numeric(gsub(pattern = '^Cluster-', replacement = '', sample.ft$cluster.id)))), sep = ''))


p <- ggplot(sample.ft,
       aes(x = cluster.id,
           y = count,
           fill = disease)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  xlab('Cluster ID') +
  ylab('Count') +
  ggtitle("Disease Profile") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  guides(fill = guide_legend(title = "Disease", ncol = 1)) +
  scale_fill_manual(values = setNames(colors.manual[1:length(unique(col.annotation$Disease))], unique(col.annotation$Disease)))
ggsave(p, 
       filename = "~/FsmPlots/DiseaseClusterEnrichmentSupport20.pdf", 
       width = 12, 
       height = 10)


## Per Cluster Plot ##
for(clust in unique(sample.ft$cluster.id)){
  local.data <- subset(sample.ft, sample.ft$cluster.id == clust)
  p <- ggplot(local.data, 
              aes(x = disease, 
                  y = count, 
                  fill = disease)) +
    geom_bar(stat = "identity", position = "dodge_2") +
    theme_minimal() +
    ggtitle(sprintf('Phenotype Enrichment (%s)', 
                    clust)) +
    theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1),
          legend.text = element_blank(),
          plot.title = element_text(hjust = 0.5, face = 'bold', size = 16)) +
    guides(fill = FALSE) + 
    xlab('') + 
    ylab('')
  ggsave(p, filename = sprintf('~/pancancer_fsm/plots/PhenotypeEnrichment%s.pdf', clust), width = 12, height = 8)
}

## Plot Network ##
custom.pathways <- read.table('~/pancancer_fsm/data/PathwayGenes.tsv',
                              header = TRUE,
                              check.names = FALSE,
                              sep = '\t',
                              stringsAsFactors = FALSE,
                              fill = TRUE)
ppi <- readRDS('~/pancancer_fsm/data/stringdb.hsa.hgnc.medconf.adj.rds')

paths <- sapply(1:ncol(custom.pathways), simplify = FALSE, function(x){
  local.genes <- na.omit(as.character(custom.pathways[,x]))
  local.genes <- local.genes[local.genes != '']
  local.genes <- local.genes[local.genes %in% colnames(ppi)]
  return(local.genes)
})
names(paths) <- colnames(custom.pathways)
rm(custom.pathways)

## Filer Matching Genes ##
all.genes <- union(fsm.edge.ft$gene.a,
                   fsm.edge.ft$gene.b)

all.genes.match <- sapply(all.genes, function(gene.id){
  local.match <- names(paths)[sapply(paths, function(x){gene.id %in% x})]
  if(length(local.match) > 0)
    return(as.character(local.match)[1])
  return('Novel')
})
all.genes.match.ft <- names(all.genes.match)[all.genes.match != 'Novel']

fsm.edge.ft$Type <- ifelse(fsm.edge.ft$gene.a %in% all.genes.match.ft | 
                             fsm.edge.ft$gene.b %in% all.genes.match.ft,
                           1, 0)
fsm.edge.ft.filt <- subset(fsm.edge.ft, fsm.edge.ft$Type == 1)
rownames(fsm.edge.ft.filt) <- 1:nrow(fsm.edge.ft.filt)

## Create Network ##
require(network)
require(sna)
require(ggplot2)
require(ggnet)

graph <- network(fsm.edge.ft.filt, 
                 directed = FALSE, 
                 hyper = FALSE, 
                 loops = FALSE, 
                 bipartite = FALSE, 
                 matrix.type = 'edgelist',
                 ignore.eval = FALSE)

graph %v% 'Pathway' <- all.genes.match[match(graph %v% 'vertex.names', 
                                             table = names(all.genes.match))]
graph %v% 'Size' <- ifelse(graph %v% 'Pathway' == 'Novel', 0.5, 1)
set.edge.attribute(graph, 'Color', ifelse(graph %e% 'Type' == 1, 'orange', 'grey40'))
set.edge.attribute(graph, 'Line', ifelse(graph %e% 'Type' == 1, 1, 2))
set.edge.attribute(graph, 'Size', ifelse(graph %e% 'Type' == 1, 0.05, 0.01))

## Plot ##
p <- ggnet2(graph,
            label = TRUE,
            label.size = 0.2,
            color = 'Pathway',
            size = 'Size',
            max_size = 1,
            alpha = 0.8,
            palette = 'Set3',
            legend.position = 'bottom',
            legend.size = 12,
            edge.size = 'Size',
            edge.lty = 'Line',
            edge.color = 'Color') +
  ggtitle('PANCANCER ATLAS FSG Network') +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold')) +
  guides(size = FALSE)
ggsave(p, filename = '~/FsmPlots/tcga.pancancer.fsg.network.support20.pdf', 
       width = 12, 
       height = 10, 
       dpi = 'retina')


## Check GBM & Thryoid Cancer Patients ##
sig.genes <- graph.cluster.list$`Cluster-63`

local.idx <- colnames(ppi) %in% sig.genes
local.ppi <- ppi[local.idx, local.idx]
local.idx <- colSums(local.ppi) == 0
local.ppi <- local.ppi[!local.idx, !local.idx]

local.graph <- network(as.matrix(local.ppi),
                       directed = FALSE)

p <- ggnet2(local.graph,
            label = TRUE,
            label.size = 2,
            color = 'orange',
            size = 6,
            edge.lty = 1,
            edge.size = 0.05,
            alpha = 0.8) +
  guides(size = FALSE) +
  ggtitle("FSG Network (Cluster-63)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 16),
        axis.title = element_blank())
ggsave(p, filename = '~/FsmPlots/Cluster63_Network.pdf', width = 12, height = 10)
