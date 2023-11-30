# Pollen Quantitative Genetics
Taking a look at the genomic variance of tomato pollen experiencing heat stress. Largely focused on implementing [kmers-based-GWAS](https://github.com/voichek/kmersGWAS/tree/master) in Nextflow.

## Contents

### kGWASflow
Testing out a Snakemake implementation of kmers-GWAS called [kGWASflow](https://github.com/akcorut/kGWASflow/tree/main).

### Nextflow
A Nextflow implementation of [kmers-based-GWAS](https://github.com/voichek/kmersGWAS/tree/master). Influenced by [kGWASflow](https://github.com/akcorut/kGWASflow/tree/main), particularily in the Conda environment setups. The key difference (other than being written in Nextflow instead of Snakemake) is that this implementation allows for some added versatility in input sequencing data. Samples can have as many libraries as you want, and you can mix long and short read libraries together. 

Some features that are coming soon:
* Ability to add different KMC parameters for short and long reads
* kmers table summary statistics
* Downstream kmer mapping and visualization, inspired by [this paper](https://doi.org/10.1002/tpg2.20374).

### R
This directory includes some GWAS data visualizations as well as a script that makes kGWASflow inputs.
