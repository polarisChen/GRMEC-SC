# GRMEC-SC
The source code and input data of GRMEC-SC

## Usage
#### Clone or download this repo.
```
git clone https://github.com/polarisChen/GRMEC-SC.git
```

#### Code structure
- ``constructSv.R```: constructing the affinity matrix S
- ```clustering_update.R```: defines the iterative steps of the model update and clustering evaluation
- ```main.py```: run the model


#### Example command
Take the datasets "cellMix" as an example.
Running environment：R version 4.2.1


## Data availability
|  Dataset              | Protocol   | Source |
| --------------------------- | ----------------------- | ----------------------- |
| ***Public***             | ***10x Multiome***      | ***https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3*** |
| ***Inhouse***          | ***CITE-seq***      | ***https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148665***     |
| ***cellMix***              | ***SNARE-seq***           | ***https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126074*** |
| ***10x-pbmc-3k***             | ***10x Multiome*** | ***https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-3-k-1-standard-2-0-0***        |
| ***Ma-2020***             | ***SHARE-seq*** | ***https://scglue.readthedocs.io/en/latest/data.html***        |
