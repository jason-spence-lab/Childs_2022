'''
BASIC SINGLE CELL ANALYSIS SCRIPT
Created by Josh Wu

Analysis used for Childs et al., 2022
Figure 3 (And Supplemental Figure 3) 

Gene list generated from Yu, Kilik, and Holloway et al., 2021, "Charting human development using a multi-endodermal organ atlas and organoid models"
DOI: 10.1016/j.cell.2021.04.028
Top 50 most differentially expressed stomach epithelial and small intestine epithelial genes used for Cell Scoring Analysis (Figure 3D-E) developed by Josh Wu
See SI_DEG.txt (small intestine genes) and stomach_DEG.txt for gene lists 

'''

from sca_run import *
#from tools.pipelines import *

figdir = #path to where you would like figures to be saved to
an_run = sca_run()
#############################################################################
## Change this to point toward your mount location for our MiStorage share ##
#############################################################################
an_run.storage_mount_point = #path to where the .h5 files are saved

## IDs of samples as represented in the metadata table
an_run.sample_list = ['299-1','299-4'] 


## List of interesting genes
an_run.add_gene_list(markers=['CDX2', 'SOX2','PDX1','TFF1','TFF2'],
					 label='Figure3B')

## Parameters used to filter the data - Mainly used to get rid of bad cells
an_run.set_filter_params(min_cells = 0, # Filter out cells 
						 min_genes = 1200, # Filter out cells with fewer genes to remove dead cells
						 max_genes = 7500, # Filter out cells with more genes to remove most doublets
						 max_counts = 50000, # Filter out cells with more UMIs to catch a few remaining doublets
						 max_mito = 0.1) # Filter out cells with high mitochondrial gene content

## Parameters used for initial clustering analysis
an_run.set_analysis_params(n_neighbors = 15, # Size of the local neighborhood used for manifold approximation
						   n_pcs = 12, # Number of principle components to use in construction of neighborhood graph
						   spread = 1, # In combination with min_dist determines how clumped embedded points are
						   min_dist = 0.4, # Minimum distance between points on the umap graph
						   resolution = 0.5,# High resolution attempts to increases # of clusters identified 
                           cell_score_lists=['SI_DEG','stomach_DEG'])#Cell Scoring Analysis for Figure 3D-E, see .txt files for gene lists

an_run.set_plot_params(size = 5, umap_obs = ['louvain','sampleName'], exp_grouping = ['louvain','sampleName'], final_quality = True) 


## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
an_run.pipe_basic(figdir)

## If you find some interesting clusters that you want to "zoom in" on and recluster, you can use the following code

# New analysis parameters for the subset of parameters
# analysis_params_ext = dict(n_neighbors = 9,
# 						n_pcs = 9,
# 						spread = 1,
# 						min_dist = 0.4,
# 						resolution = 0.4)
# an_run.pipe_ext(analysis_params_ext, figdir=figdir, extracted=['0','1','5','6','2'], load_save='adata_save.p')

