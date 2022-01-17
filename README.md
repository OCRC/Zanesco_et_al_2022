# Zanesco _et al_, 2022.
This repository contains open-source code for reproducing the analysis of single-cell RNA sequencing data exposed in the manuscript "CREB is involved in the regulation of POMC processing enzymes andα-MSH production in the hypothalamus of mice" by Zanesco _et al_, submitted to Frontiers in Endocrinology, section Neuroendocrine Science.

In this study, we harnessed the data from [Campbell et al](https://doi.org/10.1038/nn.4495), consisting of ~20,000 single-cell transcriptomes from the arcuate nucleus and the median eminence (Arc-ME), and used it to identify Creb transcriptional targets.

# Retrieve data
  We retrieved gene expression and meta-data matrices from [SingleCellPortal](https://singlecell.broadinstitute.org/single_cell/study/SCP97/a-molecular-census-of-arcuate-hypothalamus-and-median-eminence-cell-types#study-summary). You'll need to sign in with a Google account to continue:
  - [Gene expression data](https://singlecell.broadinstitute.org/single_cell/data/public/SCP97/a-molecular-census-of-arcuate-hypothalamus-and-median-eminence-cell-types?filename=expression.txt.gz)
  - [Metadata](https://singlecell.broadinstitute.org/single_cell/data/public/SCP97/a-molecular-census-of-arcuate-hypothalamus-and-median-eminence-cell-types?filename=meta.txt)

  Alternatively, download directly from the UNIX (Linux, MacOS) terminal:
  ```
  # Gene expression data
  wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93374/suppl/GSE93374_Merged_all_020816_DGE.txt.gz -O expression.txt.gz
  
  # Metadata
  wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93374/suppl/GSE93374_cell_metadata.txt.gz -O meta.txt.gz
  ```
  
  Please note that both these files must be in the same directory as the scripts (our working directory, or where you cloned the repo to), and should be named 'expression.txt.gz' and 'meta.txt.gz'.

# Single-cell analysis with dbMAP

  For this step, you'll need to have (dbMAP)[https://github.com/davisidarta/dbMAP] installed, along with (Seurat)[https://satijalab.org/seurat] for single-cell analysis in R.
  Main findings can be reproduced by executing the following line in the terminal:
  
  ```
  > R Carraro_et_al_2021_Campbell_Data_Analysis.R
  ```
  Alternatively, you may open the `Zanesco_et_al_2022_Campbell_data_analysis.R` script in your favorite IDE (i.g. RStudio) and interactively navigate the code.
  
  
# Running Arboreto and PySCENIC

  For this step, you should download some databases: https://pyscenic.readthedocs.io/en/latest/installation.html#auxiliary-datasets
  
  To run arboreto, execute the following snippets in the terminal: 
  
  
```
arboreto_with_multiprocessing.py \
    CampDataForSCENIC.tsv \
    mm_mgi_tfs.txt \
    --method grnboost2 \
    --output adj.tsv \
    --num_workers 12 \
    --seed 777
```
    
    
 And
 
 ```
pyscenic ctx adj.tsv --annotations_fname SCENIC/resources/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname CampDataForSCENIC.loom --output reg.csv --num_workers 12 --min_genes 10  databases/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather databases/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather databases/mm9-tss-centered-10kb-10species.mc9nr.feather databases/mm9-tss-centered-5kb-10species.mc9nr.feather databases/mm9-500bp-upstream-10species.mc9nr.feather
```

Additional code that could not be easily pipelined is also available at this repository.


Please send any questions and/or bug reports to davisidarta[at]fcm.unicamp.br


# Citation

TBD
