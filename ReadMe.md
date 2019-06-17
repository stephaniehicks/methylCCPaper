This repository contains the additional material and scripts for the manuscript: 

**Technology-independent estimation of cell type composition using differentially methylated regions** authored by 

* Stephanie C. Hicks, shicks19@jhu.edu
* Rafael A. Irizarry, rafa@jimmy.harvard.edu

## Data
 
All data used in this manuscript is publicly available. Here are links on how to download or access each data set: 

* [FlowSorted.Blood.450k](https://bioconductor.org/packages/release/data/experiment/html/FlowSorted.Blood.450k.html). R/Bioconductor data package containing Illumina HumanMethylation data on sorted blood cell populations. 
* [GSE95163](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95163). This dataset from [Carmona et al. (2017)](https://www.nature.com/articles/s41525-017-0012-9) measured whole blood samples on two platform technologies: (1) Illumina 450k microarray platform and (2) reduced representation bisulfite sequencing (RRBS) platform. The whole blood samples were derived from ten male individuals resulting in _N_ = 10 microarray samples and _N_=10 RRBS samples. 
* [GSE42861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42861). This dataset (_N_ = 689 samples) from [Liu et al. (2013)](https://www.ncbi.nlm.nih.gov/pubmed/23334450) studied methylation differences between rheumatoid arthritis patients and normal controls. Here, we only considered the normal controls. 
* [GSE77716](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77716). This dataset is used in [Rahmani et al. (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27018579). This is the Gala II dataset, describing a pediatric latino population with _N_ = 573 whole blood samples processed on the Illumina 450K microarray platform. There are also _N_ = 78 whole blood samples had their cell composition independently estimated using flow cytometry (processed 4 months later), which can be considered as a "gold standard".
* [BLUEPRINT Epigenome Database](http://www.blueprint-epigenome.eu). We downloaded _N_ = 44 samples from seven purified whole blood cell types, specifically B-cells, CD4T-cells, CD8T-cells, neutrophils, eosinophils, monocytes, and natural killer cells. For a given WGBS sample (e.g. CD8T-cells), we assumed the "gold standard" cell composition to be 100% CD8T-cells and 0% for the other cell types.

## methylCC software package

We used the [methylCC](https://github.com/stephaniehicks/methylCC) package to estimate the cell composition of the samples. 
