'''
BASIC SINGLE CELL ANALYSIS SCRIPT
Created by Josh Wu

Analysis used for Childs et al., 2022
Figure 2 - EGF Organoid Sample

'''

from sca_run import *
#from tools.pipelines import *

figdir = ##path to where you would like figures to be saved to
an_run = sca_run()
#############################################################################
## Change this to point toward your mount location for our MiStorage share ##
#############################################################################
an_run.storage_mount_point = #path to where the .h5 files are saved

## IDs of samples as represented in the metadata table
an_run.sample_list = ['299-4']


## List of interesting genes
an_run.add_gene_list(markers=['LGR5','OLFM4','MKI67','TOP2A','FABP2','ALPI','RBP2','BEST4','SPIB','MUC2','SPDEF','DLL1','TRPM5','TAS1R3','CHGA','NEUROD1','PAX6','ARX'],
					 label='Figure2I')


# 					 label='publicationCSC_epithelial_list')
## Parameters used to filter the data - Mainly used to get rid of bad cells
an_run.set_filter_params(min_cells = 0, # Filter out cells 
						 min_genes = 1250, # Filter out cells with fewer genes to remove dead cells
						 max_genes = 8500, # Filter out cells with more genes to remove most doublets
						 max_counts = 65000, # Filter out cells with more UMIs to catch a few remaining doublets
						 max_mito = 0.1) # Filter out cells with high mitochondrial gene content

## Parameters used for initial clustering analysis
an_run.set_analysis_params(n_neighbors = 15, # Size of the local neighborhood used for manifold approximation
						   n_pcs = 10, # Number of principle components to use in construction of neighborhood graph
						   spread = 1, # In combination with min_dist determines how clumped embedded points are
						   min_dist = 0.4, # Minimum distance between points on the umap graph
						   resolution = 0.4) # High resolution attempts to increases # of clusters identified

an_run.set_plot_params(size = 30, final_quality = True)


## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
an_run.pipe_basic(figdir)

## If you find some interesting clusters that you want to "zoom in" on and recluster, you can use the following code

# New analysis parameters for the subset of parameters
# analysis_params_ext = dict(n_neighbors = 9,
# 						n_pcs = 7,
# 						spread = 1,
# 						min_dist = 0.4,
# 						resolution = 0.3)

# an_run.pipe_ext(analysis_params_ext, figdir=figdir, extracted=['1'], load_save='adata_save.p')


