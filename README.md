# Childs_2022
This repository contains code for analyses as seen in "A Development-Inspired Niche for Homeostatic Human Mini-Intestines" - Charlie J. Childs et. al (2022).

Manuscript can be found here: https://www.biorxiv.org/content/10.1101/2022.06.12.495827v2

## Abstract
Epithelial organoids derived from intestinal tissue, also referred to as mini-intestines or mini-guts, recapitulate many aspects of the organ in vitro and can be used for biological discovery, personalized medicine, and drug development. Murine intestinal organoids represent a homeostatic system that balances stem cell maintenance within a crypt-like compartment and differentiation within a villus-like compartment. However, this homeostatic balance and spatial organization has not been achieved with human intestinal organoids. Here, we leverage single cell RNA-seq data (scRNA-seq) and high-resolution imaging to interrogate the developing human intestinal stem cell niche. We identified an EGF-family member, EPIREGULIN (EREG), as uniquely expressed in the developing crypt, and found that EREG can take the place of EGF as an in vitro niche factor. Unlike EGF, which leads to growth of thin-walled cystic organoids, EREG-organoids are spatially resolved into budded and proliferative crypt domains and a differentiated villus-like central lumen. Transcriptomics and epigenomics showed that EREG-organoids are globally similar to the native intestine while EGF-organoids have an altered chromatin landscape, downregulate the master intestinal transcription factor CDX2 and ectopically express stomach genes.

## About
Within each figure folder, the code used to create visualization found in that figure can be found. Analysis was done in either Python or R-studio mostly using the single cell RNA sequencing pipelines Scanpy (python) or Seurat (R). Major supporting packages include CellChat, scoreHIO, and Signac. 

>**Figure 1 + Supplement (Python and R)** - Single cell analysis of primary human intestine using Scanpy pipeline and Ligand-Receptor analysis using CellChat  
>**Figure 2 (Python)** - Single cell analysis of EGF-grown and EREG-grown organoids individually  
>**Figure 3 + Supplement (Python and R)** - Single cell analysis combining conditions described in Figure 2, cell scoring analysis comparing these samples to either stomach or small intestine (with accompanying gene lists used in .txt files), and reference mapping to the endoderm atlas using the scoreHIO R package developed by the Camp Lab  
>**Figure 4 + Supplement (R)** - Multiomics analysis of dual single nuclei RNA/single nuclei ATAC data of primary tissue, EREG-grown, and EGF-grown organoids.

## Scanpy Pipeline Created By: 
Joshua Wu wujos@med.umich.edu  
Mike Czerwinski czerwmj@med.umich.edu

## For Questions About Code Found in this Repo
Charlie Childs cjchilds@umich.edu

## Acknowledgments
Jason Spence Lab for support from our colleagues  
Gray Camp Lab 
Theis Lab and Satija Lab for their work on modern single cell analysis techniques
