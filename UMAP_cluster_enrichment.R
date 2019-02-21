setwd('~/pancancer_fsm/')

ppi <- readRDS('data/inbiomap.hsa.symbol.rds')
univ <- colnames(ppi)
rm(ppi)
gc()

cluster.genes <- readRDS('data/pancancer-gdc_inbiomap_umap_cluster_comparison_genes.rds')
reactome.data <- read.table('data/Ensembl2Reactome_All_Levels.txt',
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
reactome.data$V1 <- mapping$hgnc_symbol[match(reactome.data$V1, 
                                              table = mapping$ensembl_gene_id)]
reactome.list <- split(reactome.data$V1, reactome.data$V2)
reactome.list <- lapply(reactome.list, function(x){unique(na.omit(x))})
reactome.list <- reactome.list[sapply(reactome.list, length) >= 8]
names(reactome.list) <- reactome.data$V4[match(names(reactome.list),
                                               table=reactome.data$V2)]
rm(reactome.data,ensembl,mapping)
gc()

require(parallel)
cl <- makeCluster(4)
clusterExport(cl, varlist = c('univ', 'reactome.list'))
enrichRes <- parSapply(cl, 
                       cluster.genes, 
                       simplify = FALSE,
                       function(x){
                         local.enr <- sapply(reactome.list,
                                             simplify = TRUE,
                                             function(t){
                                               test.mat <- matrix(c(sum(x %in% t),
                                                                    sum(!(x %in% t)),
                                                                    sum(univ %in% t),
                                                                    sum(!(univ %in% t))),
                                                                  byrow = FALSE, ncol = 2)
                                               return(fisher.test(test.mat, alternative = 'greater')$p.value)
                                             })
                         return(p.adjust(local.enr, method = 'bonferroni'))
                       })
stopCluster(cl)
enrichRes.ft <- do.call('cbind', lapply(enrichRes, function(x){return(x)}))
write.table(enrichRes.ft, file = 'data/pancancer-gdc_inbiomap_support8_fsg_umap_cluster_enrichment.tsv',col.names = TRUE,row.names = TRUE,sep='\t',append=FALSE,quote=FALSE)
