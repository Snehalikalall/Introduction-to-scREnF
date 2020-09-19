
An Entropy Based Feature Selection- Application on Single cell RNA Sequence Data

# Summary


A novel and Robust Entropy based Feature (gene) Selection method for large Single-cell RNA-seq data, using  Renyi ´and Tsallis entropy. The SC-REnF is also demonstrated for identifying marker genes from different cell types. Our results shed some light on the single-cell clustering problem with the application of entropy based feature selection, and therefore, it will be an important tool to complement existing methods in the scRNA-seq analysis pipeline. Results are shown here in brief.


# How to use sc-REnF

## Data Loading


Load the Libraries

```
library(SingleCellExperiment)
library(edgeR)
library(scDatasets)
library(biomaRt)
library('Linnorm')
```

We will give the method demonstration with single cell RNA sequencing on 466 cells to capture the cellular complexity of the adult and fetal human brain at a whole transcriptome level. For more details about the study, see [A survey of human brain transcriptome diversity at the single cell level](https://www.pnas.org/content/112/23/7285#:~:text=Our%20results%20show%20that%20MHCI,as%20endothelial%20cells%20and%20microglia.)

Read the gene expression data using (SingleCellExperiment object), calculate CPM values and extract metadata.


```
rawdata <- readRDS("darmanis.rds")
data <- assay(rawdata)
```
For demonstration purposes, we apply a standard *Linnorm* normalization with minimum read count =5 in 10 percent cell. However any other normalization approach may be used.
Gene should be in row, Cells should be in coloumn


```
preprocessedata= normalized_data(data)

```

```
dim(preprocessedata)
[1] 10564    90

Normalized_Biase_data[1:2,1:3]
         zygote   zygote   zygote   
RPS11   8.542289 8.275977 8.199267 
ELMO2   3.419379 3.918050 3.667263 
```



