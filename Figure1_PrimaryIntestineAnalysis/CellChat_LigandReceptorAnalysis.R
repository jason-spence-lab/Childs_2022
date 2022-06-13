
# Ligand Receptor Analysis Using CellChat ---------------------------------
#Primary Intestine ligand receptor analysis for Figure 1 and Supplemental Figure 1
#Childs et al 2022

#Inital data import and download
#Set working directory and download required packages
setwd()
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

#Import single cell data from wherever it is stored
duo1271.data <- Read10X(data.dir = "/path/to/files")
duo1272.data <- Read10X(data.dir = "/path/to/files")
duo132.data <- Read10X(data.dir = "/path/to/files")
il132.data <- Read10X(data.dir = "/path/to/files")

# Initialize the Seurat object with the raw (non-normalized data).
duo1 <- CreateSeuratObject(counts = duo1271.data, project = "127DayDuo1")
duo2 <- CreateSeuratObject(counts = duo1272.data, project = "127DayDuo2")
duo3 <- CreateSeuratObject(counts = duo132.data, project = "132DayDuo")
il1 <- CreateSeuratObject(counts = il132.data, project = "132DayIl")

#Combine Seurat Objects into one merged object
combo.fetal <- merge(duo1, y = c(duo2, duo3, il1), add.cell.ids = c("127Duo1", "127Duo2", "132Duo", "132Il"), project = "fetaltissue")
combo.fetal[["percent.mt"]] <- PercentageFeatureSet(combo.fetal, pattern = "^MT-")

#Create a voilin plot to visualize raw data, filter the data, and visualize filtered data
VlnPlot(combo.fetal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
combo.fetal.clean <- subset(combo.fetal, subset = nFeature_RNA > 800 & nFeature_RNA < 3500 & percent.mt < 10 & nCount_RNA < 12000)
VlnPlot(combo.fetal.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Use SCTransform to normalize, scale, and find variable features in the dataset
combo.fetal.clean <- SCTransform(combo.fetal.clean, vars.to.regress = "percent.mt", verbose = FALSE)

#Perform dimensional reduction (PCA and UMAP)
combo.fetal.clean <- RunPCA(combo.fetal.clean, verbose = FALSE)
combo.fetal.clean <- RunUMAP(combo.fetal.clean, dims = 1:30, verbose = FALSE)
combo.fetal.clean <- FindNeighbors(combo.fetal.clean, dims = 1:30, verbose = FALSE)
combo.fetal.clean <- FindClusters(combo.fetal.clean, resolution = 0.6, verbose = FALSE)

#Visualize the data and make sure all the cell types we expect are present 
DimPlot(combo.fetal.clean, label = TRUE) + NoLegend()
DotPlot(combo.fetal.clean, features = epi_genes)
subepi <-c('LGR5','OLFM4','F3','DLL1','PDGFRA','NPY')
basic_list <-c('S100B','PLP1','STMN2','ELAVL4','CDH5','KDR','ECSCR','CLDN5','COL1A1','COL1A2','DCN','ACTA2','TAGLN','ACTG2','MYLK','EPCAM','CDH1','CDX2','CLDN4','PTPRC','HLA-DRA','ARHGDIB','CORO1A')
FeaturePlot(object = combo.fetal.clean, features = basic_list, pt.size = .1,   reduction = "umap", order = TRUE) & scale_colour_gradientn(colors=c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
FeaturePlot(object = combo.fetal.clean, features = subepi, pt.size = .1,   reduction = "umap", order = TRUE) & scale_colour_gradientn(colors=c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
FeaturePlot(object = combo.fetal.clean, features = epi_genes, pt.size = .1,   reduction = "umap", order = TRUE) & scale_colour_gradientn(colors=c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))

#Save the filtered data to your working directory
saveRDS(combo.fetal.clean, file = "cleanedfetaldata.rds")
combo.fetal.clean <- readRDS("cleanedfetaldata.rds", refhook = NULL)

#Assign cell type identities to each cluster in the UMAP (change cluster number to name of cell type)
combo.fetal.clean <- RenameIdents(combo.fetal.clean, 
                                  '0' = 'Mesenchyme',
                                  '1' = 'Mesenchyme',
                                  '2' = 'Mesenchyme', 
                                  '3' = 'Smooth Muscle',
                                  '4' = 'Mesenchyme',
                                  '5' = 'Immune',
                                  '6' = 'Immune',
                                  '7' = 'Mesenchyme',
                                  '8' = 'Subepithelial Cells',
                                  '9' = 'Mesenchyme',
                                  '10' = 'Enterocytes',
                                  '11' = 'Immune',
                                  '12' = 'Immune',
                                  '13' = 'Smooth Muscle', 
                                  '14' = 'Endothelial',
                                  '15' = 'Nerve',
                                  '16' = 'Mesenchyme',
                                  '17' = 'Immune',
                                  '18' = 'Mesenchyme',
                                  '19' = 'Nerve',
                                  '20' = 'Smooth Muscle',
                                  '21' = 'Stem Cells',
                                  '22' = 'Immune',
                                  '23' = 'Enteroendocrine Cells ',
                                  '24' = 'Mesenchyme',
                                  '25' = 'Endothelial',
                                  '26' = 'Tuft Cells',
                                  '27' = 'Goblet Cells',
                                  '28' = 'BEST4+ Enterocytes')

combo.fetal.clean$celltype <- Idents(combo.fetal.clean)

#To subset this dataset to feature only the stem cells of the epithelium (as seen in Figure 1B) 
combo.fetal.clean = subset(combo.fetal.clean, ident = c('Mesenchyme', 'Stem Cells', 'Immune', 'Endothelial', 'Smooth Muscle', 'Nerve', 'Subepithelial Cells'))
DimPlot(combo.fetal.clean)



# CellChat Analysis -------------------------------------------------------
#Follows the CellChat Vignette found here: 
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

cellchat <- createCellChat(object = combo.fetal.clean, group.by = "ident", assay = "RNA")
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat,type = "truncatedMean", trim = 0.1)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#change to pathway of interest 
pathways.show <- c("EGF") 
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
