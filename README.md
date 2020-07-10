# ExaminingPhylogenomics
Program to help examine properties of genes and supermatrices

The pupose of this is to help in examining phylogenomic matrices and understanding why an inferred relationship or difference in signal appears. The program as it stands is designed to do the investigations from the paper: [Disentangling biological and analytical factors that give rise to outlier genes in phylogenomic matrices](https://www.biorxiv.org/content/10.1101/2020.04.20.049999v1.abstract)

There are essentially four analyses that can be performed as of now. The program will perform the analyses that a given set of parameters allow. There are minimal dependencies and it is written in base python to make it easy to use. The programs that it relies on for some of the analyses are:

- [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [mafft](https://mafft.cbrc.jp/alignment/software/)
- [raxml](https://github.com/stamatak/standard-RAxML)

**If you use a function that relies on these, please, please cite the software!!!**

As of now this is the different functions the program can do. Although I am planning on updating this periodically depending on if people ask for something or I need some type of utility for investigating phylogenomic data.

To download the scripts run this command. The example data comes from this [paper](https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-10-65) where they did a great job depositing data which is incredibly helpful.

Download command

```git clone https://github.com/jfwalker/ExaminingPhylogenomics.git```

### Reverse concatenate

This is exactly what it says, it takes in a supermatrix, a partition file and splits the partitions into a folder.

To run this the following parameters must be specified:

```-d The name of an output directory```

```-z The name of the supermatrix file```

```-q The name of the RAxML formatted partition file```

**Example**

If you want you can just specify all parameters and the program will run everything at once, but for the walkthrough it's easier to do them one at a time.

After you have downloaded the folder, you'll first want to make an output folder with:

```mkdir TestFolder/```

Run the reverse concatenate with:

```python src/ExaminePhylogenomicData.py -d TestFolder/ -z example_data/ExampleConcat.fa -q example_data/ExampleGenes.model```

The folder should now be filled with the fasta's that made up the partition. To the screen will be printed what the program thought the first gene (partition) was (this is also printed to the logfile). If the first gene does not match what the partition file says, then something went wrong and you should investigate this. If the partition file cannot be broken into genes it is likely because the partition file has some type of issue or is not RAxML formatted, the logfile should have the predicted reason why it didn't work to help troubleshoot but there are so many possibilities it's hard to account for all of them.

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

**Example**

To run this you will first either need a genome file or a well sampled transcriptome file. Since this is a dataset of vertebrates and in the paper we use chicken as an example this will be used here. First you want to down load the genome for chicken from a website like [Ensembl](https://www.ensembl.org/info/data/ftp/index.html). By clicking on the link that says Fasta next to Gallus gallus you should be brought to the [ftp download](ftp://ftp.ensembl.org/pub/release-100/fasta/gallus_gallus/cds/). You can then download on the command line using wget and the link address.

```wget ftp://ftp.ensembl.org/pub/release-100/fasta/gallus_gallus/cds/Gallus_gallus.GRCg6a.cds.all.fa.gz```

Uncompress the file with

```gunzip Gallus_gallus.GRCg6a.cds.all.fa.gz```

Then you can run the orthology detection with the following command. The ```-i``` is used since Gallus is in the analysis. Also, be sure to change the path and names of the programs to the ones you have.

```python src/ExaminePhylogenomicData.py -f TestFolder/ -b blastn -m mafft -r raxmlHPC-AVX -g Gallus_gallus.GRCg6a.cds.all.fa -i Gallus```

Since Gallus is a predicted ingroup, in this case it actually is an ingroup, you should look in the IngroupOrthologyAnalysis.csv folder. In there the gene ENSGALG00000008314 that we discuss in the [paper](https://www.biorxiv.org/content/10.1101/2020.04.20.049999v1.abstract) should be flagged as True. This means to investigate further. If you open the tree, with something like [figtree](http://tree.bio.ed.ac.uk/software/figtree/) you should see a similar pattern to that which is shown as supplementary figure 5. The gene was flagged because no matter how you root it, it is impossible to have the (genomes+Gallus) defined by just two clades. Part of why I like this example is that it is imperfect, it does not detect two of the genes pointed out by [Brown and Thomson](https://academic.oup.com/sysbio/article/66/4/517/2950896), this is in part because the dataset goes back to when the ancestor of humans was on land. Blast does not identify all the genes because the similarity used is 1e-3, to identify more of the misidentified orthology you can use the protein coding sequences or re-run with a different genome. Gene duplication and loss etc. can prevent some of these issues from being detected with just one genome, so the more you check with the more you are likely to identify.


### Correcting a GWLL value to an average SSLL value

This is done because one reason the difference in likelihood can increase is based on gene length. By looking at just the gene value this only tells you the difference in contribution of genes to the final likelihood, but by looking at the Average SSLL value this can inform you that if you correct for length there is still something going on.

Running the program creates the following files:

**RAxML_perSiteLLs.Topologies_SSLL and RAxML_info.Topologies_SSLL**: This is the output from RAxML doing the SSLL topology test, where it calculates all the site likelihoods.

**GWLL.csv and AvgSSLL.csv**: This is the Gene wise log likelihood value and the average sitewise loglikelihood values, for every input species trees.

The parameters required are:

```-s This is the Species Trees```

```-z This is the supermatrix```

```-q This is the RAxML formatted partition file```

```-r This is the location of RAxML```


### Comparing the results of an ML analysis to the GWLL/SSLL analysis

The reasons these may differ is in the paper, can explain some disparity in gene tree conflict analysis

What this does is look for the different relationship among your provided species tree. If more than one exists this will proceed using the random first one it finds and the results of this should not be compared. Make sure the input species trees only have one relationship. This requires ML trees whose names can be matched to the GWLL results. The matching is identified by looking at different subsets between a period, so if your ML trees do not have a period in the name it will not be able to match. The reason a period is used is because most of these files end in .tre or some variation or they may be like RAxML and have the name after the period. To compare whether the ML results match the GWLL results you will need to have run or also be running this at the same time and give the folder of trees using the parameter:

```-t Folder of trees```

The output will be:

**Comp_tree0_tree1.csv**: This is the comparison of your first tree in the species tree file to your second tree. If you have more trees in the species tree file then it will keep doing comparisons. So three trees will give the additional files: Comp_tree0_tree2.csv and Comp_tree1_tree2.csv.

**UniqueBipartitions.txt**: How each species tree compares to eachother, the should only have one difference.

The format of the output is csv with: The gene name, Supports first tree (True/False), supports second tree (True/False), ML supports first tree (True/False/Uninformative), ML supports second tree (True/False/Uninformative).

The reason the ML may be uninformative is it possibly does not contain a bipartition that is being compared (either doesn't have sampling or tree structure). The GWLL is always forced to choose between the two topologies so it will always be True or False. 


### Finally

The program will create a logfile that is good to look over because if anything said in it seems off (It found more genes than you expected etc.) then that's something to investigate. With no specification is will create something called logfile.log or logfile can be specified with:

```-l Logfile name```

End of the day this is academic software, it attempts to automate some tedious procedures (e.g. examining every gene tree), but if something seems very off definitely contact. Also, if you have any questions or would like anything added to the software, feel free to contact for that too. You can find my email, [here](https://www.slcu.cam.ac.uk/people/walker-joseph) underneath the unfortunate picture from my first day at Cambridge.



