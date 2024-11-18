# Bacterial Mutation Rates Project
This repository is my first genomics project, created as a portfolio piece and a learning exercise in the field.

I decided to pursue an original research question, not only because it is more interesting for me, but also because
would showcase my ability to be creative and think independently. 

The long-form writeup of the project is documented in `project.pdf`, which contains a precise statement of the hypothesis, derivations and the results.

In short, the research question is:

> If one considers the genes that are essential to the survival of an organism, one would
expect that those genes are those that would least tolerate any mutations. Thus, just by looking at the mutation rates of different
genes across a population, it should be possible to infer a-priori which genes are essential for the survival of the organism.

We test this hypothesis on a sample of Escherichia Coli genomes



## Datasets
The genome data was taken from the NCBI genome datasets. We use all complete annotated e-coli genomes annotated with RefSeq. 
To download an exact copy of the dataset used to run this code, you can use the NCBIs
`datasets` command line tool, which you can install [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).

Then run the following command:
```bash
datasets download genome taxon "Escherichia coli" --include genome,cds,gtf,gff3 --annotated --assembly-level complete --assembly-source RefSeq
```

The chosen reference sequence is GCF_000005845.2 which is one of the reccomended reference
genomes on the NCBI site. 

The csv file of essential genes can be downloaded [here](https://shigen.nig.ac.jp/ecoli/pec/top.GenesListSearchAction.do?classId=1)

## Scripts

The main workflow consists of two scripts:
* `process.py`
* `analyze.py`

process.py generates a dataframe which contains mutation frequencies and other statistics for each gene. 
This process takes around 10 minutes when run on the full dataset. It produces a pickle file which is read by analyze.py

It takes a path to the downloaded NCBI dataset as well as a path to a reference genome from NCBI, which is a .fna file:

```bash

python process.py --dir=path-to-data-dir -r=path-to-reference-fna-file/cds_from_genomic.fna

```

analyze.py then reads the saved data frame and computes the various summary statistics and creates plots etc. 

The other scripts are mostly used for data inspection and debugging, and aren't necessary for the workflow:
* `process_phylo_tree.py`: We construct a phylogenetic tree so that we can remove some samples that are overly-related, probably due to duplicate uploads etc. We saved the output into duplicate_samples.txt and committed it to the repo, so that it doesn't need to be run again. Running takes ~8 hours
* `analyze_phylo_tree.py`: Takes the outputted phylo tree and finds the duplicates
* `view_data.py` is just for getting an overview of the data
* `view_gene.py` is for inspecting single genes that warrant debugging.


