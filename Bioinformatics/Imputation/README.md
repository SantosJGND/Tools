## Imputation repository.

This sub-directory holds work on imputation. As usual, this research is done with population genetic data in mind. 

The data consists a single genotype data set. Variables are variant count features ranging between 0 and 2; Samples are designed to derive from a semi-consistent population network. _Semi-consistent_ is used here to indicate that certain observations have variable pdfs, and the characteristics of the structure vary (cluster distance may change).

### Data generation

VCF files are generated using the [Genome Simulator](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Simulate_genomes/Genomic%20structure%20Simulator.ipynb) tool of the first Tools repository [link](https://github.com/SantosJGND/Tools_and_toys).

- replicated here for the specific data sets used [notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_II/blob/master/Imputation/prepare_vcfs.ipynb).


### I. Distances / Dimensionality reduction. 

Window based analysis constructs data sets of distance data with which to predict position of missing observation in incomplete data set.

> [notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_II/blob/master/Imputation/Impute_I_distances.ipynb)


### II. PCA inverse transformation.

An aside on the accuracy of PCA inverse transformation.

> [notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_II/blob/master/Imputation/INV_transform.ipynb)

### III. Cluster search. 

Dimensionality reduction and maximum likelihood cluster classification. Use for stats, imputation. 

> [notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_II/blob/master/Imputation/Reconstruct_trees.ipynb)

### Application to rice data.

**i. cluster distance and imputation**

Application of the cluster search and inference pipeline on 3000 Rice Genomes data. Focus on Japonica and cBasmati variation. Distance inference now performed within 1MB of focal target. 

> [notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_II/blob/master/Imputation/rice_ClusterSearch.ipynb)

**ii. targeted Ne estimation at local windows**

> [notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_II/blob/master/Imputation/Theta_ClusterSearch.ipynb)


![image](figures/impute_2rings2.png)