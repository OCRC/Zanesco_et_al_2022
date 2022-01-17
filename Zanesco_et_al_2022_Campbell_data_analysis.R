# Data analysis from Campbell et al. dataset, as in Zanesco et al., 2022 

# Load libraries and setwd

reticulate::use_python('/usr/bin/python3')
library(loomR)
library(reticulate)
library(Seurat)
library(SeuratWrappers)
library(plotly)

setwd("~/Documents/Bioinfo/Zanesco/")

# Download expression.txt.gz and meta.txt.gz 
# from https://singlecell.broadinstitute.org/single_cell/study/SCP97/a-molecular-census-of-arcuate-hypothalamus-and-median-eminence-cell-types


######################################################################
# Load data, create Seurat object, QC plots
######################################################################

counts <- as.matrix(read.table('expression.txt.gz', sep = '\t', header = T, row.names = 1))
meta <- read.table('meta.txt.gz', sep = '\t', header = T, row.names = 1)
dat <- CreateSeuratObject(counts = counts, meta.data = meta)

# If downloaded from GEO : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93374
# Need to subset to the filtered data shown in their paper
barcodes <- read.table('barcodes.txt')
dat <- subset(dat, cells = barcodes)

mito.genes <- grep(pattern = "^mt-", x = rownames(dat@assays$RNA@counts), value = TRUE)
percent.mito <- Matrix::colSums(dat@assays$RNA@counts[mito.genes, ])/Matrix::colSums(dat@assays$RNA@counts)
dat <- AddMetaData(object = dat, metadata = percent.mito, col.name = "percent.mito")

VlnPlot(object = dat, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"))


######################################################################
# Default SCTransform workflow; PCA and UMAP for control
######################################################################

dat <- SCTransform(dat, return.only.var.genes = F, variable.features.n = 5000)

dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)

dat <- RunUMAP(dat, dims = 1:50)
UMAPPlot(dat, group.by = 'All.Cell.Subclusters')


################################################################################ Run the diffusion basis for dbMAP
################################################################################ Import python libraries
library(reticulate)
np <- reticulate::import("numpy")
pd <- reticulate::import("pandas")
sp <- reticulate::import("scipy")
dm <- reticulate::import('dbmap')

# Deal with the matrix
data <- t(as.matrix(dat@assays$SCT@data[VariableFeatures(dat),]))
a <- r_to_py(data)
b <- sp$sparse$csr_matrix(a)

# Run the diffusion algorithm
diff = dm$diffusion$Diffusor(n_components = as.integer(200),
                             n_neighbors = as.integer(30),
                             ann_dist = as.character('cosine'),
                             n_jobs = as.integer(10),
                             kernel_use = 'simple',
                             transitions = 'True',
                             norm = 'False')
diff = diff$fit(b)
mms = diff$transform(b)
res = diff$return_dict()

# The diffusion components are the eigencomponents of the structure-weighted diffusion procedure
sc = py_to_r(res$EigenVectors)
ev = py_to_r(res$EigenValues)

################# This is to deal with the new names
rownames(sc) <- colnames(dat)
new_names <- list()
for(i in 1:length(sc)){
  new_names[i] <- paste('DC_' , as.integer(colnames(sc[i])) + 1, sep = '')
}
colnames(sc) <- as.vector(new_names)
names(ev) <- as.vector(new_names)

################# Plot the dimensionality of the dataset (similar to PCA's elbow plot)
plot(ev)

# Going back to the Seurat object
dat[["db"]] <- CreateDimReducObject(embeddings = as.matrix(sc), key = "DC_", assay = DefaultAssay(dat))

# Run UMAP for layout
um <- reticulate::import('umap')
umapper <- um$UMAP(min_dist = as.numeric(0.4), spread = as.numeric(1.5), n_epochs = as.integer(1200))
layout <- umapper$fit_transform(sc)
rownames(layout) <- colnames(dat)
plot(layout)

# Going back to the Seurat object
dat[["dbmap"]] <- CreateDimReducObject(embeddings = as.matrix(layout), key = "dbMAP_", assay = DefaultAssay(dat))
DimPlot(dat, reduction = 'dbmap', group.by = 'All.Cell.Clusters', pt.size = 0.5)

###############################################################################
# Producing the plots
###############################################################################

# Creb1 and Pomc
FeaturePlot(dat, reduction='dbmap', 
            features = c('Creb1', 'Pomc'), pt.size = 1, order = T)

# Pomc region 
FeaturePlot(dat, reduction='dbmap', 
            features = c('Ppp1r9a', 'Npy5r', 'Nptx2', 'Aak1', 'Chl1'), pt.size = 1, order = T)

# Major Creb1 transcriptional targets in the Arc-ME
FeaturePlot(dat, reduction='dbmap', 
            features = c('Ism1', 'Atp6v1c2', 'Ccnd2', 'Ccdc39', 'Stc1', 'Mgat3', 'Gmcl1', 'Pik3r2', 'Utp20', 'Asic1', 'Usp45', 'Hspbap1'), pt.size = 1, order = T)




##############################################################################
# SCENIC
##############################################################################
# Write neurons expression values for SCENIC analysis of transcriptional networks
write.table(as.matrix(dat@assays$RNA@counts), file = 'CampDataForSCENIC.tsv', row.names = T, col.names = T)

# Also save as Loom for downstream SCENIC analysis
dat_loom <- create('CampDataForSCENIC.loom', data=dat@assays$RNA@counts)

# In terminal
# GRNs (arboreto):
# arboreto_with_multiprocessing.py     CampDataForSCENIC.loom    mm_mgi_tfs.txt     --method grnboost2     --output adj.tsv     --num_workers 12     --seed 777

# CTX:
# pyscenic ctx adj.tsv --annotations_fname /home/davi/Documents/Bioinfo/SCENIC/resources/motifs-v9-nr.mgi-m0.001-o0.0.tbl
# --expression_mtx_fname CampDataForSCENIC.loom --output reg.csv --num_workers 12
# /home/davi/Documents/Bioinfo/SCENIC/databases/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
# /home/davi/Documents/Bioinfo/SCENIC/databases/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
# /home/davi/Documents/Bioinfo/SCENIC/databases/mm9-tss-centered-10kb-10species.mc9nr.feather
# /home/davi/Documents/Bioinfo/SCENIC/databases/mm9-tss-centered-5kb-10species.mc9nr.feather
# /home/davi/Documents/Bioinfo/SCENIC/databases/mm9-500bp-upstream-10species.mc9nr.feather


# AUCell:
# pyscenic aucell CampDataForSCENIC.loom reg.csv --num_workers 11 -o scenic_auc.loom


# Reading the output

adj <- read.table('adj.tsv', header = T)

adj_creb <- adj[adj$TF == 'Creb1',]
adj_creb_filtered <- adj_creb[adj_creb$importance > 1,]

write.table(adj_creb_filtered, file = 'adj_creb_filtered.tsv', row.names = T, col.names = T)


# load the library
library(forcats)
library(ggplot2)
library(dplyr)

# Reorder following the value of another column:
name = adj_creb_filtered$target
val = adj_creb_filtered$importance
adj_creb_filtered %>%
  mutate(name = fct_reorder(name, val)) %>%
  ggplot( aes(x=name, y=val)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()










