library(umap)

setwd('~/pancancer_fsm/')

## FSG Matrix
fsg.mat <- readRDS('data/pancancer-gdc_fsg_dysregulation_mat.rds')
fsg.mat <- scale(t(fsg.mat))

## UMAP ##
custom.config = umap.defaults
custom.config$n_components <- 2
custom.config$n_neighbors <- 20
custom.config$verbose <- TRUE
custom.config$n_epochs <- 1000

umap.res <- umap(fsg.mat, method = 'naive', config = custom.config)

## Save
saveRDS(umap.res, file = 'data/pancancer-gdc_fsg_dysregulation_mat_umap.rds')
