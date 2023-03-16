
### [fastq management](fastq_management/)

Fastq filters.

### [contig bridging using reads](bridge_contigs/)

Map reads to assembly, find reads that bridge contigs.Find *contig connected components*; summary statistics.

### [Uniform distribution tests](Distribution/)

Kolmogorov-Smirnof test across contigs, readme inside. 

### [gff distances](gff_distances/)

Distances and overlap between segments in separate gff files. README inside./ 

### [orthologous gene identification](GeneXtract_v1/)

Map and fasta to list of assemblies, extract best positions and corresponding sequences, 
- subset orfs; 
- predict protein;
- construct graph; 

### [Inversion Informative Markers](IIMs/)

Calculate ld within inversions, use Inversion matrix as input. 


### [Merging Annotations](merge_PAVs/)

Merge presence absence calls from different software. Given standardized input format. 

i. merge Insertions and Deletions separately.
ii. Produce overlap statistics and 3-way venn diagrams.

### [Demographic inference](Ocoarcata_demography/)

Demographic inference using single nucleotide polimorphism data and fastsimcoal2.6 software.

- Three population scenarios, with and without migration;
- Setup for independent deployment across subgenomes.

Used for inference of Oryza Coarctata.


### [Inversion typing](INVtype/)

This was meant to determine inversion orientatio from mapping strands. But it was a bust. Still.

Given region in reference assembly (samtools format), and mapping bamfile (against reference assembly):
- extract alignments to selecte region.
- re-map against extract fasta to generate paf.

- readme.md inside. 


## Factor Process

Generalizing data analysis.

>  [notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_II/blob/master/Factor_process/Factor_walk.ipynb)

## Factor Dict

Recursive hierarchic index dictionary creation from factor array.

> [script](Factor_dict/factor_dict.py)

## Imputation

_VCF, dimensionality reduction, distances, KDE_

> [directory](https://github.com/SantosJGND/Tools_II/tree/master/Imputation)







