'''
Alignment Outlying Signal
'''
import sys
import Extras
import FastaTools
import argparse
import LikelihoodTests
import TreeTools

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
    parser.add_argument("-i", "--iqtree_location", required=False, type=str, help="""Location of iqtree""")
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
        dna = FastaTools.blast_it(args.FolderOfFastas,args.blast,args.genome_file,outw)
        #2nd) Extract the sequences and create new ones with genome seqs in it
        FastaTools.extract_blast("Genome.fa","seq_ortho_test.rawblast","TempOrthoFolder/")
        #3rd) Align the new fastas
        FastaTools.align_folder("TempOrthoFolder/",args.LocationOfMafft,outw)
        #4th) Infer Tree from the new alignment
        LikelihoodTests.RunRaxml("TempOrthoFolder/",args.raxml_location,dna,outw)
        #5th) Examine placement of genome compared to rest
        TreeTools.identify_ortho_issue("TempOrthoFolder/",outw)
    
    '''
    If nothing ran then print this
    '''
    if NoRun == True:
        Extras.get_time("Program Died :(", outw)
        print "To run the program needs the files and locations of programs necessary to perform the analyses needed. Each setting will say in it's description what it is used for.\n\nFor info run: python OutlierDissection.py -h\n"

if __name__ == "__main__":
    main()


