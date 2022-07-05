# Basic Multiomics Analysis Script ----------------------------------------
#This analysis leverages Seurat and Signac pipelines to analyze dual single nuclei RNA/ATAC multiomic data
#Samples in manuscript include fresh crypts, EREG and EGF grown human organoids
#This script was created by Charlie Childs 2022
#Childs et al 2022 Figure 4 and Supplemental Figure 4
#Raw files created from the 10x Cell Ranger ARC analysis pipeline from the Michigan Advanced Genomics Core


# Import Multiomic Data and Create One Object -----------------------------
#This section imports all data as individual objects and then merges these objects into one
#Preprocessing for all parts of figure 4 and supplemental figure 4
#Set working directory
setwd('/set/path/to/directory')

#Download Packages
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse) #not necessary but trying to write files
library(RColorBrewer)
library(leaflet)
library(GenomicRanges)
library(future)
set.seed(1234)

plan('multiprocess', workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

#Import peak .bed files for each sample (file = ata_peaks.bed)
peaks.fc <- read.table(
  file = '/path/to/atac_peaks.bed',
  col.names = c('chr', 'start', 'end')
)

peaks.ereg <- read.table(
  file = '/path/to/atac_peaks.bed/',
  col.names = c('chr', 'start', 'end')
)

peaks.egf <- read.table(
  file = '/path/to/atac_peaks.bed/',
  col.names = c('chr', 'start', 'end')
)

#Create a common set of peaks between all objects for seamless merging into one object
gr.fc <- makeGRangesFromDataFrame(peaks.fc)
gr.ereg <- makeGRangesFromDataFrame(peaks.ereg)
gr.egf <- makeGRangesFromDataFrame(peaks.egf)
combined.peaks <- IRanges::reduce(x=c(gr.fc,gr.egf,gr.ereg))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

#Read in cell barcodes for each sample
md.fc <- read.table(
  file = '/path/to/per_barcode_metrics.csv',
  stringsAsFactors = FALSE,
  sep = ',',
  header = TRUE,
  row.names = 1
)[-1, ]

md.ereg <- read.table(
  file = '/path/to/per_barcode_metrics.csv',
  stringsAsFactors = FALSE,
  sep = ',',
  header = TRUE,
  row.names = 1
)[-1, ]

md.egf <- read.table(
  file = '/path/to/per_barcode_metrics.csv',
  stringsAsFactors = FALSE,
  sep = ',',
  header = TRUE,
  row.names = 1
)[-1, ]

#Initial filtering of low count cells
md.fc <- md.fc[md.fc$passed_filters > 300, ]
md.ereg <- md.ereg[md.ereg$passed_filters > 300, ]
md.egf <- md.egf[md.egf$passed_filters > 300, ]

#Create fragment objects for each sample (file = )
frag_fc <- CreateFragmentObject('/path/to/atac_fragments.tsv.gz')
frag_ereg <- CreateFragmentObject('/path/to/atac_fragments.tsv.gz')
frag_egf <- CreateFragmentObject('/path/to/atac_fragments.tsv.gz')

#Create matrix for each sample of peaks x cells
fc_counts <- FeatureMatrix(
  fragments = frag_fc,
  features = combined.peaks,
  cells = rownames(md.fc)
)

ereg_counts <- FeatureMatrix(
  fragments = frag_ereg,
  features = combined.peaks,
  cells = rownames(md.ereg)
)

egf_counts <- FeatureMatrix(
  fragments = frag_egf,
  features = combined.peaks,
  cells = rownames(md.egf)
)

#Genome annotation for h38 
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- 'UCSC'

#Create a Seurat Object for each sample with proper counts and frag path
#Fresh crypts RNA object
fc_counts <- Read10X_h5('/volumes/umms-spencejr/01_RAW_RNASEQ_AGC_SHARE/4757-CC/10x_analysis_4757-CC/Sample_4757-CC-1/filtered_feature_bc_matrix.h5')
fc_matrix <- CreateSeuratObject(
  counts = fc_counts$`Gene Expression`,
  assay = 'RNA'
)
#Add Fresh Crypts ATAC object
fc_matrix[['ATAC']] <- CreateChromatinAssay(
  counts = fc_counts$Peaks,
  sep = c(':', '-'),
  fragments = frag_fc,
  annotation = annotation
)

#1ng/ml EREG sample RNA object
ereg_counts <- Read10X_h5('/volumes/umms-spencejr/01_RAW_RNASEQ_AGC_SHARE/4757-CC/10x_analysis_4757-CC/Sample_4757-CC-2/filtered_feature_bc_matrix.h5')
ereg_matrix <- CreateSeuratObject(
  counts = ereg_counts$`Gene Expression`,
  assay = 'RNA'
)
#Add 1ng/ml EREG sample ATAC object
ereg_matrix[['ATAC']] <- CreateChromatinAssay(
  counts = ereg_counts$Peaks,
  sep = c(':', '-'),
  fragments = frag_ereg,
  annotation = annotation
)

#100ng/ml EGF sample RNA object
egf_counts <- Read10X_h5('/volumes/umms-spencejr/01_RAW_RNASEQ_AGC_SHARE/4757-CC/10x_analysis_4757-CC/Sample_4757-CC-4/filtered_feature_bc_matrix.h5')
egf_matrix <- CreateSeuratObject(
  counts = egf_counts$`Gene Expression`,
  assay = 'RNA'
)
#Add 100ng/ml EGF sample ATAC object
egf_matrix[['ATAC']] <- CreateChromatinAssay(
  counts = egf_counts$Peaks,
  sep = c(':', '-'),
  fragments = frag_egf,
  annotation = annotation
)

#Give each sample an identity before merging so cells from these samples can be tracked
egf_matrix$dataset <- '100ng/ml EGF'
ereg_matrix$dataset <- '1ng/ml EREG'
fc_matrix$dataset <- 'Fresh Crypts'

#Merge all Seurat Objects into one
combo_matrix <- merge(
  x = fc_matrix,
  y = list(ereg_matrix, egf_matrix),
  add.cell.ids = c('Fresh Crypts', '1ng/ml EREG', '100ng/ml EGF')
)

#Save merged object 
saveRDS(combo_matrix, file = 'combo_matrix.rds')
combo_matrix

#Use this line to switch back and forth between the RNA and ATAC aspects of the object throughout the script
DefaultAssay(combo_matrix) <- 'ATAC' 

#Compute per-cell quality control metrics from ATAC data 
combo_matrix <- NucleosomeSignal(combo_matrix)
combo_matrix <- TSSEnrichment(combo_matrix)

#Pre-filter violin plot visualization 
VlnPlot(
  object = combo_matrix,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  group.by = 'dataset',
  ncol = 4,
  pt.size = 0
)

#Filtering step 
combo_matrix <- subset(
  x = combo_matrix,
  subset = nCount_ATAC < 45000 &
    nCount_RNA < 15000 &
    nCount_ATAC > 100 &
    nCount_RNA > 100 &
    TSS.enrichment > .5
)

#Post processing voilin plot
VlnPlot(
  object = combo_matrix,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  group.by = 'dataset',
  ncol = 4,
  pt.size = 0
)

#Save filtered data object
saveRDS(combo_matrix, file = 'combo_matrix')


# Peak Calling ------------------------------------------------------------
#Call peaks using MACS2
peaks <- CallPeaks(combo_matrix, macs2.path = 'path/to/macs2')
peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
saveRDS(peaks, file = 'peaks.rds')

#Quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(combo_matrix), #combined fragment object from all fragment files
  features = peaks, #combined peaks
  cells = colnames(combo_matrix) #unfiltered cells do I need filtered?
)
#Save Macs2counts
saveRDS(macs2_counts, file = 'macs2_counts.rds')

#Create a new assay using the MACS2 peak set and add it to the Seurat object
combo_matrix[['peaks']] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = list(frag_fc,frag_ereg,frag_egf), #tried list of fragpaths, Fragments(filtered...),
  annotation = annotation
)

saveRDS(combo_matrix, file = 'filtered_combo_matrix.rds')
combo_matrix <- readRDS('filtered_combo_matrix.rds', refhook = NULL)


# Data Visualization ------------------------------------------------------
#RNA analysis, normalization, dimensional reduction, and UMAP visualization
DefaultAssay(combo_matrix) <- 'RNA'
combo_matrix <- SCTransform(combo_matrix, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

#ATAC analysis
DefaultAssay(combo_matrix) <- 'ATAC'
combo_matrix <- RunTFIDF(combo_matrix)
combo_matrix <- FindTopFeatures(combo_matrix, min.cutoff = 'q0')
combo_matrix <- RunSVD(combo_matrix)
combo_matrix <- RunUMAP(combo_matrix, reduction = 'lsi', dims = 2:50, reduction.name = 'umap.atac', reduction.key = 'atacUMAP_')

#Calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities
combo_matrix <- FindMultiModalNeighbors(combo_matrix, reduction.list = list('pca', 'lsi'), dims.list = list(1:50, 2:50))
combo_matrix <- RunUMAP(combo_matrix, nn.name = 'weighted.nn', reduction.name = 'wnn.umap', reduction.key = 'wnnUMAP_')
combo_matrix <- FindClusters(combo_matrix, graph.name = 'wsnn', resolution = 0.2, algorithm = 3, verbose = FALSE)

#Plot RNA UMAP, ATAC UMAP, and WNN UMAP (multiomic data), Supplemental Figure 4A-B
p1 <- DimPlot(combo_matrix, reduction = 'umap.rna', group.by = 'dataset', cols = c('#6D6E71', '#00AEEF', '#991B45'), repel = TRUE) + ggtitle('RNA')
p2 <- DimPlot(combo_matrix, reduction = 'umap.atac', group.by = 'dataset', cols = c('#6D6E71', '#00AEEF', '#991B45'), repel = TRUE) + ggtitle('ATAC')
p3 <- DimPlot(combo_matrix, reduction = 'wnn.umap', group.by = 'dataset', cols = c('#6D6E71', '#00AEEF', '#991B45'), repel = TRUE) + ggtitle('WNN')
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

#Fresh Crypt sample will have some contaminating mesenchyme, so filter for epithelial cells (CDH1+, VIM -)
extracted_epi <- subset(combo_matrix, sct_CDH1 > 0 & sct_VIM==0)
#Visualize epithelial only dataset as we did with the entire dataset, Figure 4A
p1 <- DimPlot(extracted_epi, reduction = 'umap.rna', group.by = 'dataset', cols = c('#6D6E71', '#00AEEF','#991B45'), repel = TRUE) + ggtitle('RNA')
p2 <- DimPlot(extracted_epi, reduction = 'umap.atac', group.by = 'dataset', cols = c('#6D6E71', '#00AEEF','#991B45'), repel = TRUE) + ggtitle('ATAC')
p3 <- DimPlot(extracted_epi, reduction = 'wnn.umap', group.by = 'dataset', cols = c('#6D6E71', '#00AEEF','#991B45'), repel = TRUE) + ggtitle('WNN')
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

# Motif Analysis Using Signac ---------------------------------------------
#Figure 4B and 4C
#Packages required for this section
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(presto)

#Genome used for mapping the data to has the scaffolds named differently to the BSgenome. This chunk of code fixs this issue.
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(extracted_epi))) %in% main.chroms)
extracted_epi[['ATAC']] <- subset(extracted_epi[['ATAC']], features = rownames(extracted_epi[['ATAC']])[keep.peaks])

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object 
DefaultAssay(extracted_epi) <- 'ATAC'
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(extracted_epi), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
extracted_epi <- SetAssayData(extracted_epi, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

#Compute a per cell motif activity score 
extracted_epi <- RunChromVAR(
  object = extracted_epi,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

markers_rna <- presto:::wilcoxauc.Seurat(X = extracted_epi, group_by = 'dataset', assay = 'data', seurat_assay = 'SCT')
markers_motifs <- presto:::wilcoxauc.Seurat(X = extracted_epi, group_by = 'dataset', assay = 'data', seurat_assay = 'chromvar')
motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0('RNA.', colnames(markers_rna))
colnames(markers_motifs) <- paste0('motif.', colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(extracted_epi, id = motif.names)

#Visualize both gene expression and motif score for each cell for genes/associated motifs of interest
motif.name <- ConvertMotifID(extracted_epi, name = c('CDX2','GATA4', 'HNF4A'))
gene_plot <- DotPlot(extracted_epi, assay = 'SCT', features = c('CDX2','GATA4', 'HNF4A'), group.by = 'dataset')& scale_colour_gradientn(colors=c('lightgrey','#FFFFD9','#EDF8B1','#C7E9B4','#7FCDBB','#41B6C4','#1D91C0','#225EA8','#253494','#081D58'))
motif_plot <- DotPlot(extracted_epi, assay = 'chromvar', features = motif.name, group.by = 'dataset')& scale_colour_gradientn(colors=c('lightgrey','#FFFFD9','#EDF8B1','#C7E9B4','#7FCDBB','#41B6C4','#1D91C0','#225EA8','#253494','#081D58'))
gene_plot | motif_plot

# Link Peaks and Genes ----------------------------------------------------
#Figure 4D and Supplemental Figure 4C and 4D

DefaultAssay(extracted_epi) <- 'peaks'
extracted_epi <- RegionStats(extracted_epi, genome = BSgenome.Hsapiens.UCSC.hg38)
extracted_epi <- addGCBias(extracted_epi, genome = BSgenome.Hsapiens.UCSC.hg38)

# link combo_peaks to genes
extracted_epi <- LinkPeaks(
  object = extracted_epi,
  peak.assay = 'peaks',
  expression.assay = 'SCT',
  genes.use = c('CDX2','GATA4', 'HNF4A')
)

#Visualize linked gene peak plots
CoveragePlot(
  object = extracted_epi,
  region = 'CDX2',
  features = 'CDX2',
  expression.assay = 'SCT',
  group.by =  'dataset',
  extend.upstream = 1000,
  extend.downstream = 1000
)