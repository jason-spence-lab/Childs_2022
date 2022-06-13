# Endoderm Atlas Reference Mapping (scoreHIO R Package) -------------------
#This script contains analysis code for Figure 3F-G of Childs et al 2022
#It leverages the scoreHIO package developed by the Camp Lab
#Original paper: Yu, Kilik, Holloway et al, Cell, 2021 
#DOI: 10.1016/j.cell.2021.04.028

#Intial data import and filtering
# Load the package
library(devtools)
devtools::install_github('Camp-Lab/scoreHIO')
library(scoreHIO)
library(Seurat)
library(ggplot2)
library(reticulate)
library(uwot)
library(RANN)
library(quadprog)

#load data for each sample
oneereg.data <- Read10X(data.dir = "path/to/files")
egf.data <- Read10X(data.dir = "path/to/files")

# Initialize the Seurat Objects for each sample
oneereg <- CreateSeuratObject(counts = oneereg.data, project = "1ereg")
egf <- CreateSeuratObject(counts = egf.data, project = "egf")

#Merge the two Seurat Objects into one Object 
combo.ent <- merge(oneereg, y = egf, add.cell.ids = c("1EREG", "EGF"), project = "enteroidcombo")
combo.ent[["percent.mt"]] <- PercentageFeatureSet(combo.ent, pattern = "^MT-")

#Create a voilin plot to visualize raw data, filter the data, and visualize filtered data
VlnPlot(combo.ent, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
combo.ent.clean <- subset(combo.ent, subset = nFeature_RNA > 1200 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA < 50000)
VlnPlot(combo.ent.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Normalize the data
combo.ent.norm<- NormalizeData(combo.ent.clean)


# scoreHIO Package --------------------------------------------------------
#From here, analysis follows scoreHIO documentation which can be found here:
#https://github.com/Camp-Lab/scoreHIO
#Output represents figure 3F-G

# You could estimate the organ fidelity of these query HIOs by running this:
seurat_object <- score_fidelity(
  que_obj = combo.ent.norm,
  organ_ref_dir = "path/to/Ref_data_for_projection_to_fetal_atlas/",
  group_by = "orig.ident", 
  return_object = TRUE, 
  output_dir="path/to/output", 
  probability_file_name="Plot_barplot_projected_developing_organ_probability.pdf",
  do_projection=TRUE,
  projection_file_name="Plot_UMAP_RSS_projection_to_developing_reference.pdf")
