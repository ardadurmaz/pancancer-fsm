library(Rtsne)
library(Matrix)
library(stringi)
library(parallel)
library(umap)

setwd('~/pancancer_fsm/')

ppi <- readRDS('data/inbiomap.hsa.symbol.rds')
univ <- colnames(ppi)
result.wd <- '~/pancancer_fsm/PANCANCER-GDC/fsm_results/'

### Read FSM Results ##
cl <- makeCluster(20)
fsm.res <- parSapply(cl,
                     list.files(path = result.wd,
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
                         sample.name <- as.character(stri_match(fName, regex = 'TCGA-.{2}-.{4}-\\d{2}.{1}', mode = 'first'))
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
rownames(fsm.res) <- 1:nrow(fsm.res)
gc()

## Extract Genes per FSG ##
local.s <- stri_match_all(fsm.res$fsg, regex = "\\d\\:([^\\(]+?)\\]")
cl <- makeCluster(20)
fsm.genes.pp <- parSapply(cl, local.s, simplify = TRUE, function(x){
  require(stringi)
  return(paste(sort(unique(unlist(strsplit(x[,2],split = '.',fixed = TRUE)))), collapse = ','))
})
stopCluster(cl)
fsm.res$fsg <- fsm.genes.pp
write.table(fsm.res, file = 'data/gdc_pancancer/PANCANCER-GDC_InBioMap_Support8_FSM.tsv',
            col.names = TRUE,row.names = FALSE,sep='\t',append=FALSE,quote=FALSE)
fsm.genes.pp.uniq <- unique(fsm.genes.pp)
saveRDS(fsm.genes.pp.uniq, file = 'data/gdc_pancancer/PANCANCER-GDC_InBioMap_Support8_FSM_Unique.rds')

rm(local.s)
gc()
## Map Expression Values ##
expr <- readRDS('~/pancancer_fsm/data/gdc_pancancer/gdcExprNorm.rds')
expr <- expr[,colnames(expr) %in% unique(fsm.res$sample.id)]
expr.scaled <- t(scale(t(expr)))

fsm.mat <- matrix(0, 
                  ncol = ncol(expr.scaled), 
                  nrow = length(fsm.genes.pp.uniq), 
                  dimnames = list(paste('fsg', 1:length(fsm.genes.pp.uniq), sep = ''),
                                  colnames(expr.scaled)))

for(f in 1:nrow(fsm.mat)){
  cat(sprintf('\rProcessing FSG: %d', f))
  fsg.genes <- unlist(strsplit(fsm.genes.pp.uniq[f], split = ','))
  fsg.scores <- round(sqrt(colSums(expr.scaled[which(rownames(expr.scaled) %in% fsg.genes),]^2)), digits = 3)
  fsm.mat[f,] <- fsg.scores
}


## t-SNE ##
tsne.res <- Rtsne(scale(t(fsm.mat)), 
                  dims = 2, 
                  perplexity = 20, 
                  theta = 0.1, 
                  pca = TRUE, 
                  partial_pca = TRUE,
                  pca_center = FALSE,
                  pca_scale = FALSE,
                  normalize = FALSE,
                  num_threads = 20)


## UMAP ##
custom.config = umap.defaults
custom.config$n_components <- 2
custom.config$n_neighbors <- 20
custom.config$verbose <- TRUE
custom.config$n_epochs <- 1000

umap.res <- umap(scale(t(fsm.mat)), method = 'naive', config = custom.config)

plot.data <- data.frame('Dimension.A' = umap.res$layout[,1],
                        'Dimension.B' = umap.res$layout[,2],
                        'Disease' = pheno.data$disease_type[match(rownames(umap.res$data),
                                                                  table=pheno.data$sample)],
                        stringsAsFactors = FALSE)

ggplot(plot.data,
       aes(x = Dimension.A,
           y = Dimension.B,
           color = Disease)) +
  geom_point(
    shape = 20,
    size = 2) +
  theme_minimal() +
  theme(legend.direction = 'vertical',
        legend.box = 'vertical') +
  ggtitle('PANCANCER UMAP') +
  guides(color = guide_legend(nrow = 27,
                              override.aes = list(size = 4)))


## Plot ##
library(ggplot2)
pheno.data <- read.table('~/pancancer_fsm/data/gdc_pancancer/GDC-PANCAN.GDC_phenotype.tsv',
                         header = TRUE,
                         sep = '\t',
                         comment.char = '',
                         quote = '',
                         stringsAsFactors = FALSE)

plot.data <- data.frame('Dimension.A' = tsne.res$Y[,1],
                        'Dimension.B' = tsne.res$Y[,2],
                        'Disease'  = pheno.data$disease_type[match(colnames(fsm.mat),
                                                                   table=pheno.data$sample)])
annot.data.x <- tapply(plot.data$Dimension.A, plot.data$Disease, median)
annot.data.y <- tapply(plot.data$Dimension.B, plot.data$Disease, median)
annot.data.text <- names(annot.data.x)

p <- ggplot(plot.data,
            aes(x = Dimension.A,
                y = Dimension.B,
                color = Disease)) +
  geom_point(
    shape = 20,
    size = 1,
    alpha = 0.8) +
  theme_minimal() +
  theme(legend.direction = 'vertical',
        legend.box = 'vertical') +
  ggtitle('PANCANCER UMAP') +
  guides(color = FALSE) + #guide_legend(nrow = 27, override.aes = list(size = 4))) +
  annotate("text", 
           x = annot.data.x,
           y = annot.data.y,
           label = annot.data.text,
           size = 3,
           alpha = 0.6,
           fontface = 2)

write.table(plot.data, file = '~/FsmPlots/PANCANCER-GDC_InBioMap_Support8_FSG_tSNE.tsv', col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, quote = FALSE)

## Genes Only ##
all.fsm.genes <- Reduce(union,
                        sapply(fsm.genes.pp.uniq, simplify = FALSE, function(x){
                          unlist(strsplit(x,split=',',fixed=TRUE))
                        }))
expr.filt <- expr[match(all.fsm.genes, table = rownames(expr)),]
expr.filt.scaled <- scale(t(expr.filt))
gc()

tsne.res <- Rtsne(expr.filt.scaled, 
                  dims = 2, 
                  perplexity = 20, 
                  theta = 0.01, 
                  pca = FALSE, 
                  pca_center = FALSE, 
                  pca_scale = FALSE, 
                  partial_pca = FALSE, 
                  normalize = FALSE, 
                  num_threads = 20)
plot.data <- data.frame('Dimension.A' = tsne.res$Y[,1],
                        'Dimension.B' = tsne.res$Y[,2],
                        'Disease'  = pheno.data$disease_type[match(colnames(fsm.mat),
                                                                   table=pheno.data$sample)])
annot.data.x <- tapply(plot.data$Dimension.A, plot.data$Disease, median)
annot.data.y <- tapply(plot.data$Dimension.B, plot.data$Disease, median)
annot.data.text <- names(annot.data.x)


p <- ggplot(plot.data,
            aes(x = Dimension.A,
                y = Dimension.B,
                color = Disease)) +
  geom_point(
    shape = 20,
    size = 1,
    alpha = 0.8) +
  theme_minimal() +
  theme(legend.direction = 'vertical',
        legend.box = 'vertical') +
  ggtitle('PANCANCER t-SNE (Expression Only)') +
  guides(color = FALSE) + #guide_legend(nrow = 27, override.aes = list(size = 4))) +
  annotate("text", 
           x = annot.data.x,
           y = annot.data.y,
           label = annot.data.text,
           size = 3,
           alpha = 0.6)
