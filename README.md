# ExaminingPhylogenomics
Program to help examine properties of genes and supermatrices

The pupose of this is to help in examining phylogenomic matrices and understanding why an inferred relationship or difference in signal appears. The program as it stands is designed to do the investigations from the paper: [Disentangling biological and analytical factors that give rise to outlier genes in phylogenomic matrices](https://www.biorxiv.org/content/10.1101/2020.04.20.049999v1.abstract)

There are essentially four analyses that can be performed as of now. The program will perform the analyses that a given set of parameters allow. There are minimal dependencies and it is written in base python to make it easy to use. The programs that it relies on for some of the analyses are:

- [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [mafft](https://mafft.cbrc.jp/alignment/software/)
- [raxml](https://github.com/stamatak/standard-RAxML)

**If you use a function that relies on these, please, please cite the software!!!**

As of now this is the different functions the program can do. Although I am planning on updating this periodically depending on if people ask for something or I need some type of utility for investigating phylogenomic data.

### Reverse concatenate

This is exactly what it says, it takes in a supermatrix, a partition file and splits the partitions into a folder.

To run this the following parameters must be specified:

```-d The name of an output directory```

```-z The name of the supermatrix file```

```-q The name of the RAxML formatted partition file```

### Identify gene trees that are possibly misidentified orthology

This will take a genome (Or well sampled transcriptome but just called genome from now on) in fasta format, find which of the orthologs the sequences from the genome match with blast, align the orthologs with the genome sequences, create a tree and look for a patten in the tree that would indicate the gene has misidentified orthology. It creates some output files and an output folder.

The files are:

**OutgroupOrthologyAnalysis.csv**: This is the results of if there were any red flags noticed in the positions of the genomes assuming they were outgroups. The file is delimited as: GeneName (from folder specified), gene name with genomes (in TempOrthoFolder), whether a red flag appeared based on the placement of those genomes, what the largest bipartition of genome seqs was.

**IngroupOrthologyAnalysis.csv**: This is the results of if there were any red flags noticed in the positions of the genomes assuming they were ingroups. The file is delimited as: GeneName (from folder specified), gene name with genomes (in TempOrthoFolder), whether a red flag appeared based on the placement of those genomes, what the reason for the flag is or is not.

**gene_conversion.csv**: This is a comma separated file of what your gene names were in the folder vs. what they are now

**Conversion.csv**: This is what the gene names were in the input genome file vs. what they are in the Genome.fa file.

**Genome.fa**: This is a file where the names have been converted.

The folder is:

**TempOrthoFolder**: This contains the new fastas, alignments and inferred gene trees with the genome seqs included.

The reason two orthology files are created is because if the genome you use is from a species that should be outgroup to everything in the analysis then you want to look for a different patten then the one for ingroup. The pattern for outgroup would be that all the genome sequences fall into one bipartition. The pattern you look for in the ingroup case is that they cannot be defined by two bipartitions. The reason for two is that you may root on one and have the other be inside your ingroup clade, so two bipartitions are ok, any more than that are an issue. For what this pattern may show up as, the paper above has it in the supplemental data. Finally, if the genome you choose is part of the ingroup then that you should specify the name of this in the gene trees using ```-i```.

The other parameters this analysis requires are:

```-f The folder where the input seqs are```

```-b the path to blastn or blastp depending on molecule type you are analyzing```

```-m the path to mafft```

```-r the path to raxml classic```

```-g A genome or transcriptome  or some kind of well sampled fasta```


### Correcting a GWLL value to an average SSLL value

This is done because one reason the difference in likelihood can increase is based on gene length. By looking at just the gene value this only tells you the difference in contribution of genes to the final likelihood, but by looking at the Average SSLL value this can inform you that if you correct for length there is still something going on.

Running the program creates the following files:

**RAxML_perSiteLLs.Topologies_SSLL and RAxML_info.Topologies_SSLL**: This is the output from RAxML doing the SSLL topology test, where it calculates all the site likelihoods.

**GWLL.csv and AvgSSLL.csv**: This is the Gene wise log likelihood value and the average sitewise loglikelihood values, for every input species trees.

**Comp_tree0_tree1.csv** This is the comparison of your first tree in the species tree file to your second tree. If you have more trees in the species tree file then it will keep doing comparisons. So three trees will give the additional files: Comp_tree0_tree2.csv and Comp_tree1_tree2.csv.










