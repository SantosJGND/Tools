### PopGen/

Population genetics applications.


### [Population genetics jobs](PopGen/)

- Fst calculatons,
- Window-based estimates of diversity (pi and tajD),
- Window-based conversion from plink tyo;
- treemix replicates consensus tree extraction ;
- Demogaphic simulations using SLiM - alternative 2 population scenarios with and without positive a$
- ePSMC demographic inference and selfing-rate estimation ;
- ABC inference;


### FST/
- Fst calculatons,

scripts:
> fst_calc.sh; plot_launch.sh;

### PopGenStats/

- Window-based estimates of diversity (pi and tajD) from VCF file.
	- takes repeat mask file. 
	- calculations per sub-group: see IDs/ directory.
scripts:
> pi_get.sh

### ePSMC/


### treemix/

Extracting consensus tree from treemix replicates in python. 

### RealData/

Given plink file and list of subset IDs (underscore delimited); 
	- extract local genomic windows of given size, convert to MS format. 

*used to compare with SLiM simulations in SLiM/; 

### SLiM/

Demogaphic simulations using SLiM - alternative 2 population scenarios with and without positive and background selection.

subdirectories:
- `CAL/`: calculate statistics per simulation ; 
- `ePSMC`: demographic inference on simulated dta. 
- `recipes`: slim recipes: demographic scenarios and genome structure (gene density)
- `environments`: yaml files for environments used.
- `pop1_CrashNEUT`: example run directory. 

### ABC/

Approximate Bayesian Computation: determining best demographic parameters for single model for future model compaison. 

subdirectories: 
- `programs`: copy of main scripts, utilities scripts. 
- `tempsubtrop.v3` : example run. 


### PCA_Genomcics/ 

Tools for Genomic Analysis reliant on Principal Component Analysis (PCA) transformation of haplotype data, including Local Ancestry Inference. 

subdirectories: 
- 'AMOVA': AMOVA calculation through PCA, comparison with traditional methods and other metrics. 
- 'PCA_Fst': Calculate Fst from euclidean distances. 
- 'Simulate_genomes': Using PCA to simulate structure. 
- 'VCF_Analysis': Diversity, Structure, Classification, Simulation and Local Ancestry Inference.
- 'Simulate_genomes': Simulate genomes with varying loclal genomic structure (e.g. precise simulation of introgression).
- 'Cluster_shape': Calculate Allele frequency distribution, correcting for local introgression.



