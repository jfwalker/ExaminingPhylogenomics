import sys
import os
import Extras

#return an array in the order of the partition file with
#the lenght of each partition, assumes partition file was
#checked earlier on for formatting
def get_parts_array(parts):
    
    p = open(parts,"r")
    array = []
    part_array = []
    name_array = []
    
    for line in p:
        
        line = line.strip("\r\n")
        array = line.split()
        start = int(array[3].split("-")[0])
        stop = int(array[3].split("-")[1])
        #getting the complex math going on here
        gene_len = (stop - start + 1)
        part_array.append(gene_len)
        name_array.append(array[1])
    
    return name_array,part_array

#turn a fasta to hash/map/dictionary/associative array or whatever: {name => fasta}
def fasta_to_hash(supermatrix):

    file = open(supermatrix,"r")
    
    name = ""
    HASH = {}
    
    for line in file:
        line = line.strip("\r\n")
        if line[0] == ">":
            name = line
            HASH[name] = ""
        else:
            HASH[name] += line.upper()

    return HASH

#This is a sanity check of the partition file
#should be RAxML format
def get_part_format(first_line,logfile):


    array = []
    dna = True
    start = ""
    stop = ""

    first_line = first_line.strip("\n\r")
    
    if len(first_line.split(",")) != 2:
        
        print "This doesn't seem to be raxml format, missing the comma or too many commas"
        logfile.write("This doesn't seem to be raxml format, missing the comma or too many commas\n")
        print "Example: DNA, partition_name = Start-Stop"
        logfile.write("Example: DNA, partition_name = Start-Stop\n")
        Extras.get_time("Died from partition file format",logfile)
        sys.exit(0)
   
    if len(first_line.split("=")) != 2:
        
        print "This doesn't seem to be raxml format, missing the equal sign or too many equal signs"
        logfile.write("This doesn't seem to be raxml format, missing the equal sign or too many signs\n")
        print "Example: DNA, partition_name = Start-Stop"
        logfile.write("Example: DNA, partition_name = Start-Stop\n")
        Extras.get_time("Died from partition file format",logfile)
        sys.exit(0)

    if len(first_line.split(" ")) != 4:
    
        print "Something is off about the partition file but cant diagnose what it is"
        logfile.write("Something is off about the partition file but cant diagnose what it is\n")
        print "Example: DNA, partition_name = Start-Stop"
        logfile.write("Example: DNA, partition_name = Start-Stop\n")
        Extras.get_time("Died from partition file format",logfile)
    
    else:
        array = first_line.split()
        
        if array[0][:-1] == "DNA":
            dna = True
        else:
            dna = False
        
        start = array[3].split("-")[0]
        stop = array[3].split("-")[1]
        print "===================================================================="
        print "First line of the partition file says the data is DNA\nThe gene name is: " + array[1] + " \nThe start position is " + start + "\nThe stop position is " + stop + "\nIf any of this seems off, don't trust the results"
        logfile.write("====================================================================\n")
        logfile.write("First line of the partition file says the data is DNA\nThe gene name is: " + array[1] + " \nThe start position is " + start + "\nThe stop position is " + stop + "\nIf any of this seems off, don't trust the results\n")
        logfile.write("====================================================================\n")
        print "===================================================================="
        
        return dna

#Determine if the user kept the /
def check_folder(file):

    if file[-1] == "/":
        return True
    else:
        return False
    


#function for turning supermatrix into corresponding gene fastas
def rev_cat(folder,supermatrix,parts,logfile):
    
    #memory efficiency is for losers
    Supermatrix_HASH = {}
    Supermatrix_HASH = fasta_to_hash(supermatrix)
    
    file = open(parts,"r")
    #sanity check the partition file and find out if it's dna or aa
    
    has_slash = False
    
    has_slash = check_folder(folder)
    print has_slash
    
    #core of rev_cat
    array = []
    start = 0
    stop = 0
    seq = ""
    first_line = True
    
    for line in file:
    
        line = line.strip("\n\r")
        if first_line == True:
            dna = get_part_format(line,logfile)
            first_line = False
        
        #0: ignore, 1: gene name, 2: ignore, 3: location
        array = line.split()
        start = int(array[3].split("-")[0]) - 1
        stop = int(array[3].split("-")[1])
        
        if has_slash == True:
            gene_out = folder + array[1] + ".fa"
            gene_fasta = open(gene_out,"w")
        else:
            gene_out = folder + "/" + array[1] + ".fa"
            gene_fasta = open(gene_out,"w")
        
        for i in Supermatrix_HASH:
            
            missing = 0
            seq = Supermatrix_HASH[i][start:stop]

            if dna == True:
                #get missing data for dna
                missing = seq.count("N") + seq.count("-")
                if missing != (stop - start):
                    gene_fasta.write(i + "\n" + Supermatrix_HASH[i][start:stop] + "\n")
            else:
                
                missing = seq.count("X") + seq.count("-")
                if missing != stop:
                    gene_fasta.write(i + "\n" + Supermatrix_HASH[i][start:stop] + "\n")
    
    Extras.get_time("Divided into folders",logfile)

#This is just because some fastas can be super weird
def genome_format(genome,logfile):
    
    outw = open("Genome.fa","w")
    outw2 = open("Conversion.csv","w")
    outw2.write("Was,Is\n")
    g = open(genome,"r")
    count =  1
    for line in g:
        line = line.strip("\r\n")
        if line[0] == ">":
            outw.write(">Genome" + str(count) + "\n")
            outw2.write(line + "," + ">Genome" + str(count) + "\n")
            count += 1
        else:
            outw.write(line + "\n")
    Extras.get_time("New genome file is called Genome.fa",logfile)
    Extras.get_time("Table of previous name to current name is Conversion.csv",logfile)
    outw.close()
        
#Gets all the sequences from a folder and condenses them into one file
def check_and_comp_seqs(Folder,logfile):
    
    array = []
    array = os.listdir(Folder)
    has_slash = check_folder(Folder)
    if has_slash == False:
        Folder += "/"
    
    Folder2 = "TempOrthoFolder/"
    os.system("mkdir TempOrthoFolder/")
    conv = open("gene_conversion.csv","w")
    out_seq = open("seq_ortho_test.fa","w")
    count = 0
    conv.write("Was,Is\n")
    
    for x in array:
    
        count += 1
        count2 = 0
        file = open(Folder + x,"r")
        outw = open(Folder2 + "gene" + str(count) + ".fa","w")
        conv.write(x + "," + "gene" + str(count) + "\n")
        
        for line in file:
            
            #unalign the sequence
            line = line.strip("\r\n")
            
            if line[0] == ">":
                outw.write(line + "\n")
                if count2 < 10:
                    out_seq.write(">gene" + str(count) + "_" + str(count2) + "\n")
                    
            else:
                line = line.replace("-","")
                outw.write(line + "\n")
                if count2 < 10:
                    out_seq.write(line + "\n")
                count2 += 1
    
    Extras.get_time("Genes converted",logfile)
    logfile.write("New File seq_ortho_test.fa made for blast\n")
    logfile.write("New File gene_conversion.csv to keep track of old and new gene names\n")
    conv.close()
    file.close()
    outw.close()
    out_seq.close()



#blast a folder of fastas against a genome
def blast_it(Folder,blast,genome,logfile):

    #check nuc or prot
    if blast[-1] == "n":
        Extras.get_time("The data appears to be nucleotide",logfile)
        cmd = "makeblastdb -in Genome.fa -out Genome.fa -dbtype 'nucl'"
        cmd2 = "blastn -db Genome.fa -query seq_ortho_test.fa -evalue 1e-3 -num_threads 2 -max_target_seqs 100 -out seq_ortho_test.rawblast -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'"
        dna = True
    elif blast[-1] == "p":
        Extras.get_time("The data appears to be protein",logfile)
        cmd = "makeblastdb -in Genome.fa -out Genome.fa -dbtype 'prot'"
        cmd2 = "blastp -db Genome.fa -query seq_ortho_test.fa -evalue 1e-3 -num_threads 2 -max_target_seqs 100 -out seq_ortho_test.rawblast -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'"
        dna = False
    else:
        print "Can't tell what this is nuc or aa? please have blast either be blastn or blastp"
        logfile.write("Died at blast check, please have blast named blastp or blastn\n")
        sys.exit(0)
    
    #make sure the genome is formatted in an easy way to use
    #a file called genome.fa is
    print "Making genome file easier to use"
    #genome_format(genome,logfile)
    
    #Check seqs from folder for being aligned
    check_and_comp_seqs(Folder,logfile)
    
    print "Making blast database with command: " + cmd
    logfile.write("Making blast database with command " + cmd + "\n")
    
    os.system(cmd)
    
    print "Running blast database with command: " + cmd2
    logfile.write("Making blast database with command " + cmd2 + "\n")
    os.system(cmd2)
    
    return dna

#takes in a blast file and returns what is associated with what assuming multiple
#hits were found
def process_blast(blast_file):

    array = []
    HASH = {}
    
    b = open(blast_file,"r")
    for line in b:
        line = line.strip("\r\n")
        array = line.split()
        gene = array[0].split("_")[0]
        
        if gene in HASH:
            if array[2] not in HASH[gene]:
                HASH[gene].append(array[2])
        else:
            HASH[gene] = []
            HASH[gene].append(array[2])
    return HASH

#Get a list of which blast are associated with what
#Then put those into respective genes
def extract_blast(genome,blast,Folder):

    HASH = {}
    GenomeHASH = {}
    #find out which ones are associated with what
    HASH = process_blast(blast)
    
    #Figure out if there are duplicates from the Folder
    array = []
    array = os.listdir(Folder)
    
    print "Getting genome into hash"
    GenomeHASH = fasta_to_hash(genome)
    
    for x in array:
        name = x.split(".")[0]
        if name in HASH:
            #no duplicate so toss it
            if len(HASH[name]) == 1:
                os.system("rm " + Folder + x)
            #has more than one so add it
            else:
                outw = open(Folder+x,"a")
                for y in HASH[name]:
                    outw.write(">" + y + "\n" + GenomeHASH[">"+y] + "\n")
                outw.close()
        else:
            os.system("rm " + Folder + x)
        
def align_folder(Folder,mafft,logfile):

    Extras.get_time("Aligning files in folder",logfile)
    array = []
    array = os.listdir(Folder)
    
    for x in array:

        cmd = str(mafft) + " --auto --maxiterate 1000 TempOrthoFolder/" + x + " > " + "TempOrthoFolder/" + x + ".aln"
        logfile.write("Mafft command is: " + cmd)
        os.system(cmd)
    Extras.get_time("Finished aligning",logfile)
