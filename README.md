# Childs_2022
This repository contains code for analyses as seen in "EPIREGULIN creates a developmental niche for spatially organized human intestinal enteroids" - Charlie J. Childs et. al (2023) JCI Insights DOI:10.1172/jci.insight.165566.
<p align="center">
<img width="308" alt="EREGorganoid" src="https://user-images.githubusercontent.com/55200067/173452750-33f59291-4e88-43e2-9971-20c86046bcd3.png">
</p>

## Abstract
Epithelial organoids derived from intestinal tissue, also referred to as mini-intestines or mini-guts, recapitulate many aspects of the organ in vitro and can be used for biological discovery, personalized medicine, and drug development. Murine intestinal organoids represent a homeostatic system that balances stem cell maintenance within a crypt-like compartment and differentiation within a villus-like compartment. However, this homeostatic balance and spatial organization has not been achieved with human intestinal organoids. Here, we leverage single cell RNA-seq data (scRNA-seq) and high-resolution imaging to interrogate the developing human intestinal stem cell niche. We identified an EGF-family member, EPIREGULIN (EREG), as uniquely expressed in the developing crypt, and found that EREG can take the place of EGF as an in vitro niche factor. Unlike EGF, which leads to growth of thin-walled cystic organoids, EREG-organoids are spatially resolved into budded and proliferative crypt domains and a differentiated villus-like central lumen. Transcriptomics and epigenomics showed that EREG-organoids are globally similar to the native intestine while EGF-organoids have an altered chromatin landscape, downregulate the master intestinal transcription factor CDX2 and ectopically express stomach genes.

## About This Repo 
Within each figure folder, the code used to create visualization for that specific figure will be housed. Analysis was done in either Python or R-studio mostly using the single cell RNA sequencing pipelines Scanpy (python) or Seurat (R). Major supporting packages include CellChat, scoreHIO, and Signac. Parameters used for analysis in the manuscript can be found in the specific script used for that figure. 

>**Figure 1 + Supplement (Python and R)** - Single cell analysis of primary human intestine using Scanpy pipeline and Ligand-Receptor analysis using CellChat  
>**Figure 2 (Python)** - Single cell analysis of EGF-grown and EREG-grown organoids individually  
>**Figure 3 + Supplement (Python and R)** - Single cell analysis combining conditions described in Figure 2, cell scoring analysis comparing these samples to either stomach or small intestine (with accompanying gene lists used in .txt files), and reference mapping to the endoderm atlas using the scoreHIO R package developed by the Camp Lab  
>**Figure 4 + Supplement (R)** - Multiomics analysis of dual single nuclei RNA/single nuclei ATAC multiomic data of primary tissue, EREG-grown, and EGF-grown organoids. K-means clustering and motif analysis completed by the Verzi Lab at Rutgers Univeristy.  

## Data Availability
Sequencing data generated and used by this study is deposited at EMBL-EBI ArrayExpress. Datasets for human fetal intestine (ArrayExpress: E-MTAB-9489, https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9489/, and previously published work for our lab (Holloway_Czerwinski_Tsai_2020)), scRNA-seq of human fetal intestinal organoids (ArrayExpress: E-MTAB-11912, https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11912/), and dual snRNA/snATAC-seq of human fetal intestinal organoids (ArrayExpress: E-MTAB-11908, https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11908/) have been deposited. 

## Datasets, Scripts, and Analysis Created By: 
Charlie Childs cjchilds@umich.edu  
(For questions about code found in this repo, please contact Charlie)

Scanpy Pipeline Created By Joshua Wu and Mike Czerwinski 

## Acknowledgments
Jason Spence Lab   
Gray Camp Lab    
Mike Verzi Lab  
Theis Lab and Satija Lab for their work on modern single cell analysis techniques
