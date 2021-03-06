'''
Alignment Outlying Signal
'''
import sys
import Extras
import FastaTools
import argparse
import LikelihoodTests
import TreeTools
import os

def generate_argparser():

    parser = argparse.ArgumentParser(
        prog="ExaminePhylogenomicData.py",
        )
    parser.add_argument("-s", "--species_trees", required=False, type=str, help="""
    Alternative Topologies for analysis. Tests [Alignment Length, Rate Heterogeneity vs. ML tree conflict]""")
    parser.add_argument("-t", "--treeset", required=False, type=str, help="""
    Folder of gene trees [Rate Heterogeneity vs. ML tree conflict]""")
    parser.add_argument("-z", "--supermatrix", required=False, type=str, help="""
    Supermatrix for checking rate influence and/or alignment length""")
    parser.add_argument("-q", "--partition_file", required=False, type=str, help="""
    Partition file""")
    parser.add_argument("-b", "--blast", required=False, type=str, help="""Location of blastn or blastp depending on molecule type""")
    parser.add_argument("-f", "--FolderOfFastas", required=False, type=str, help="""These are used to check orthology""")
    parser.add_argument("-m", "--LocationOfMafft", required=False, type=str, help="""Used for orthology check""")
    parser.add_argument("-g", "--genome_file", required=False, type=str, help="""Used for orthology check""")
    parser.add_argument("-r", "--raxml_location", required=False, type=str, help="""Location of raxml""")
    parser.add_argument("-i", "--ingroup", required=False, type=str, help="""Genome of the ingroup for orthology""")
    parser.add_argument("-d", "--divide_matrix", required=False, type=str, help="""divide supermatrix into the genes and place them in specified folder""")
    parser.add_argument("-l", "--log_file", type=str, help="""Output file""")
    return parser


def main(arguments=None):

    parser = generate_argparser()
    args = parser.parse_args(arguments)
    
    if args.log_file:
        outw = open(args.log_file,"a")
    else:
        outw = open("logfile.log","a")

    Extras.get_time("Starting", outw)
    
    NoRun = True
    
    #reverse concatenate
    if args.divide_matrix and args.supermatrix and args.partition_file:
        Extras.get_time("Dividing supermatrix into genes",outw)
       
        NoRun = False
        FastaTools.rev_cat(args.divide_matrix,args.supermatrix,args.partition_file,outw)
    
    '''
    Alignment based outlier detection
    '''
    #performs GWLL and calculates average SSLL
    if args.species_trees and args.supermatrix and args.partition_file and  args.raxml_location:
    
        NoRun = False
        Extras.get_time("=================Alignment Test=======================",outw)
        Extras.get_time("Species trees and supermatrix identified, running two top test using RAxML",outw)
        LikelihoodTests.RunRaxML_TwoTop(args.species_trees,args.supermatrix,args.partition_file,args.raxml_location,outw)
        '''
        Heterogeneity examination (requires the GWLL so only runs post GWLL
        '''
        if args.treeset:
            Extras.get_time("=================Heterogeneity Test=======================",outw)
            #1) read in the deltaGWLL file HASH{GeneName => [test1,test2,test3,...]}
            HASH_GWLL = {}
            HASH_GWLL = Extras.get_deltaGWLL("deltaGWLL.csv")
            #2) Get the bipartitions that differ between the species trees
            Array_Tree = []
            Array_Tree = TreeTools.get_dif_biparts(args.species_trees,outw)
            #3) Compare which bipartition is in the gene set
            TreeTools.get_trees_matching_species(Array_Tree,HASH_GWLL,args.treeset,outw)
            
            

    '''
    Orthology based outlier detection
    '''
    #blasts genome of an ingroup to the genes and identifies if
    #the genome splits the ingroups into multiple bipartitions
    if args.FolderOfFastas and args.blast and args.LocationOfMafft and  args.raxml_location and args.genome_file:
        
        NoRun = False
        Extras.get_time("=================Orthology Test=======================",outw)
        #1st) blast sequences in folder against the genome
        isFile = Extras.check_if_happened("TempOrthoFolder/","TempOrthoFolder already exists moving to Genome Extract")
        if isFile == False:
            dna = FastaTools.blast_it(args.FolderOfFastas,args.blast,args.genome_file,outw)
            #2nd) Extract the sequences and create new ones with genome seqs in it
            FastaTools.extract_blast("Genome.fa","seq_ortho_test.rawblast","TempOrthoFolder/")
        #3rd) Align the new fastas
        isFile = Extras.check_if_happened_in_dir("TempOrthoFolder/","aln","end","alignment happened moving to trees")
        if isFile == False:
            FastaTools.align_folder("TempOrthoFolder/",args.LocationOfMafft,outw)
        #4th) Infer Tree from the new alignment
        isFile = Extras.check_if_happened_in_dir("TempOrthoFolder/","RAxML","beg","trees happened moving to ortho examine")
        if isFile == False:
            dna = False
            if args.blast[-1] == "n":
                dna = True
            LikelihoodTests.RunRaxml("TempOrthoFolder/",args.raxml_location,dna,outw)
        #5th) Examine placement of genome compared to rest
        isFile = Extras.check_if_happened("IngroupOrthologyAnalysis.csv","Orthology check already completed\nTo re-run the process please move or remove the following from the directory:\n-TempOrthoFolder/\n-Genome.fa\n-IngroupOrthologyAnalysis.csv\n-OutgroupOrthologyAnalysis.csv\n-seq_ortho_test.rawblast")
        if isFile == False:
            if args.ingroup:
                TreeTools.identify_ortho_issue("TempOrthoFolder/",args.ingroup,outw)
            else:
                TreeTools.identify_ortho_issue("TempOrthoFolder/","",outw)
    
    '''
    If nothing ran then print this
    '''
    if NoRun == True:
        Extras.get_time("Program Died :(", outw)
        print "To run the program needs the files and locations of programs\nnecessary to perform the analyses needed. Each setting\nwill say in it's description what it is used for.\n\nFor info run: python ExaminePhylogenomicData.py -h\n"
        print "================================================================"
        print "If you want to extract the fastas from your supermatrix and put\nthem into a separate folder you'll need to give the following arguments\n -d the folder you want the results in\n -z The supermatrix \n -q The partition file"
        print "================================================================"
        print "For the orthology test you will need to give the following arguments\n -f FolderWithFastas\n -b Location of blast\n -m Location of mafft\n -r Location of raxml\n -g Fasta of GenomeFile or reference of some kind\n -i If the genome you are using is one of your ingroups, specify which ingroup it is"
        print "================================================================"
        print "For alignment based outlying behavior, you'll need\n -s SpeciesTrees \n -z Supermatrix.fa \n -q PartitionFile (raxml formatted) \n -r Location of raxml"
        print "================================================================"
        print "If you want to compare the ML topology to the GWLL then you'll want to also specify\n -t Folder of gene trees"
        print "================================================================"
        print "If none of this makes sense or you would like anything added to\nthe program please open an issue on github or contact jfwalker@umich.edu"

if __name__ == "__main__":
    main()


