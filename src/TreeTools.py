import sys
import os
import Extras
'''
Tools for investigating the trees, lots of random old code from the old
days in here
'''

class Node:
    def __init__(self):
        self.label = ""
        self.length = 0.0
        self.time_length = 0.0
        self.parent = None
        self.children = []
        self.data = {}
        self.istip = False
        self.height = 0
        
    def add_child(self,child):
        #make sure that the child is not already in there
        assert child not in self.children
        self.children.append(child)
        child.parent = self

#This takes in the newick and the
#seq data then puts them in a data
#structure that can be preorder or
#postorder traversed pretty easily
def build(instr):
	#print "Entered build"
	root = None
	name_array =[]
	index = 0
	nextchar = instr[index]
	begining = "Yep"
	keepgoing = True
	current_node = None
	#keeps going until the value becomes false
	while keepgoing == True:
		#This situation will only happen at the very beginning but
		#when it hits this it will create a root and change begining
		#to no
		if nextchar == "(" and begining == "Yep":
				
			root = Node()
			current_node = root
			begining = "No"
		#This happens anytime their is an open bracket thats not the
		#beginning
		elif nextchar == "(" and begining == "No":
		
			newnode = Node()
			current_node.add_child(newnode)
			current_node = newnode
		#This indicates that you are in a clade and tells the 
		#program to move back one to grab the sister to the clade
		elif nextchar == ',':
		
			current_node = current_node.parent
		#This says you are closing a clade and therefore it moves
		#back to where the parent node is which allows the name
		#to be added to the parent node
		elif nextchar == ")":
			#print "Closing Clade"
			current_node = current_node.parent
			index += 1
			nextchar = instr[index]
			while True:
			
				if nextchar == ',' or nextchar == ')' or nextchar == ':' \
					or nextchar == ';' or nextchar == '[':
					break
				name += nextchar
				index += 1
				nextchar = instr[index]
			current_node.label = name
			index -= 1
		#This indicates everything is done so keepgoing becomes false
		elif nextchar == ';':
		
			keepgoing = False
			break
		#This indicates you have branch lengths so it grabs the branch
		#lengths turns them into floats and puts them in the current node
		elif nextchar == ":":
			index += 1
			nextchar = instr[index]
			while True:
				if nextchar == ',' or nextchar == ')' or nextchar == ':' \
					or nextchar == ';' or nextchar == '[':
					break
				branch += nextchar
				index += 1
				nextchar = instr[index]
			current_node.length = float(branch)
			index -= 1
		#This is for if anywhitespace exists
		elif nextchar == ' ':
		
			index += 1
			nextchar = instr[index]
		#This is for when any taxa name is hit, it will concatenate
		#the taxa names together and add the name
		else: # this is an external named node
		
			newnode = Node()
			current_node.add_child(newnode)
			current_node = newnode
			current_node.istip = True
			while True:
				if nextchar == ',' or nextchar == ')' or nextchar == ':' \
					or nextchar == ';' or nextchar == '[':
					break
				name += nextchar
				index += 1
				nextchar = instr[index]
			current_node.label = name
			name_array.append(name)
			index -= 1
		if index < len(instr) - 1:
			index += 1
		nextchar = instr[index]
		name = ""
		branch = ""
	return name_array,root

#get the other side of the bipartition
def get_right(clade_names, all_names):

    mis1 = list(set(all_names) - set(clade_names))
    clade_names.append("|")
    return ",".join(clade_names + mis1)

#get segments of a bipart
def clade_post_order(clade,clade_names):
    
    for x in clade.children:
        if x.istip:
            clade_names.append(x.label)
        clade_post_order(x,clade_names)
    return clade_names

#Post order traverse the whole tree
def post_order(tree,all_names,t_to_clade):
    
    for x in tree.children:
        #account for trees that don't have support
        if x.children and x.label == "":
            #print "Clade does not have support value"
            clade_names = []
            clade = []
            clade_names = clade_post_order(x,clade_names)
            clade = get_right(clade_names, all_names)
            t_to_clade.append(clade)
            
        elif x.children:
            #print "Clade has support value: " + x.label
            clade_names = []
            clade = []
            clade_names = clade_post_order(x,clade_names)
            clade = get_right(clade_names, all_names)
            t_to_clade.append(clade)
            
        post_order(x,all_names,t_to_clade)
    return t_to_clade

def matches_and_taxa(genome_seqs,clade):

    if(all(x in clade for x in genome_seqs)):
        return len(clade)
    else:
        return 99999999

#get the largest bipartition that is only genomes
#if that is equal to all the genomes then there is
#no outgroup that's popped into the ingroup
def get_smallest_outgroup(clades,genome_seqs):
    
    good = True
    left_array = []
    right_array = []
    largest = []
    temp = 0
    
    for x in clades:
        left_array = x.split("|")[0][:-1].split(",")
        good = only_genomes(left_array,genome_seqs)
        if good == True:
            if temp < len(left_array):
                temp = len(left_array)
                largest = left_array
        right_array = x.split("|")[1][1:].split(",")
        good = only_genomes(right_array,genome_seqs)
        if good == True:
            if temp < len(right_array):
                temp = len(right_array)
                largest = right_array
    
    if temp == len(genome_seqs):
        return largest,True
    else:
        return largest,False

#checks to see if a bipartitions is composed of only genomes
def only_genomes(bip,genomes):
    if(set(bip).issubset(set(genomes))):
        return True
    else:
        return False

def get_number_of_unique_arrays_in_array_of_arrays(further_test):
    
    #grabbed from stack overflow because everything has been asked on there
    unique_data = [list(x) for x in set(tuple(x) for x in further_test)]
    return unique_data

#takes in two arrays and removes the target sequences from
#Array1
def remove_overlap(Array1,SequencesYouWantRemoved):

    return_array = []
    for x in Array1:
        if x not in SequencesYouWantRemoved:
            return_array.append(x)
    return return_array

#get the number of one arrays seqs that are in another
def Array_matches(Array1,Array2,number_to_beat):

    count = 0
    for x in Array1:
        if x in Array2:
            count += 1
    if count == number_to_beat:
        return count
    else:
        return 0


#get clades without genomes and those with genomes
def get_smallest_ingroup_given(clades,genome_seqs,ingroupname):
    
    genome_seqs.append(ingroupname)
    largest_genome = []
    largest_genome,ans = get_smallest_outgroup(clades,genome_seqs)
    
    #If they all form one giant cluster or if they separate into one cluster
    #and one spare indicating they exhibit no more than 2 edges
    if len(largest_genome) == len(genome_seqs):
        return "Forms one clade",False
    if len(largest_genome) == (len(genome_seqs)-1):
        return "Forms one clade and one tip",False
    
    #if only two taxa exist, then it's possible one is ingroup and one is outgroup
    #hard to tell if it's any orthology issue
    if len(genome_seqs) == 2:
        return "Separate biparts but only two genomes so no orthology statement",False
    
    genomes_remaining = remove_overlap(genome_seqs,largest_genome)
    
    number_to_beat = 0
    second_largest = []
    #remove and look for next largest (do they again take up two edges)
    for x in clades:
        left_array = x.split("|")[0][:-1].split(",")
        right_array = x.split("|")[1][1:].split(",")
        new_right = remove_overlap(right_array,largest_genome)
        new_left = remove_overlap(left_array,largest_genome)
        #get the bipart with the most only genomes
        genomes_left_bipart = only_genomes(new_left,genomes_remaining)
        if genomes_left_bipart == True:
            if number_to_beat < len(new_left):
                number_to_beat = len(new_left)
                second_largest = new_left
        
        genomes_right_bipart = only_genomes(new_right,genomes_remaining)
        if genomes_right_bipart == True:
            if number_to_beat < len(new_right):
                number_to_beat = len(new_right)
                second_largest = new_right
    
    running_genomes = []
    running_genomes.append(largest_genome)
    running_genomes.append(second_largest)
    #get the size of the contents of running_genomes since it is an array of arrays
    running_genomes_true_size = len(largest_genome) + len(second_largest)
    #This means that they've formed two perfect clades, which is also
    #indistinguishable from proper orthology in an unrooted tree
    if running_genomes_true_size == len(genome_seqs):
        return "Forms two clades",False

    return "Hard to define by two edges",True

'''
#get smallest Genome seq bipartition
#This is a bit trickier if using an ingroup
#Conditions would be that the ingroup can form
#two clades, one that only contains itself duplication
#prior to orthology and one that contains the rest plus other seqs
#Their cannot be a third one though
def get_smallest_ingroup(clades,genome_seqs):
    
    left_array = []
    right_array = []
    bip_in_question = []
    further_test = []
    other_clades = []
    smallest_clade = 0
    compassing_clade = 9999999
    only_genome_clade = []
    only_genome_clade_count = 0
    unique_genome_clades = []
    genome_seq_pos = 0
    all_clades_of_genome = []
    
    
    for y in genome_seqs:
        temp =[]
        temp.append(y)
        bip_in_question,ans = get_smallest_outgroup(clades,temp)
        ans = only_genomes(bip_in_question,genome_seqs)
        if ans == True:
            only_genome_clade_count += 1
            bip_in_question.sort()
            only_genome_clade.append(bip_in_question)
        else:
            bip_in_question.sort()
            further_test.append(bip_in_question)
    other_clades = get_number_of_unique_arrays_in_array_of_arrays(further_test)
    unique_genome_clades = get_number_of_unique_arrays_in_array_of_arrays(only_genome_clade)

    #If unique_genome_clades has more than 1 that can be a red flag
    #if other clades has more than 1 that can be a red flag
    #basically if the sum of the two of those is great than 2 than thats
    #a violation of orthology on the tree
    genome_seq_pos = len(unique_genome_clades) + len(other_clades)
    all_clades_of_genome = unique_genome_clades + other_clades
    
    #current issue (Genome,(Genome,Seq)) is counted as two
    #Maybe check to see if the other clades all contain the same Seq(s)?
    #Other issue Genome over counting
    if 2 < genome_seq_pos:
        return True,all_clades_of_genome
    else:
        return False,all_clades_of_genome
'''
def get_seqs_from_genomes(array):
    
    array2 = []
    for x in array:
        if x[0:6] == "Genome":
            array2.append(x)
    return array2

#designed to see if all genomes are in one connected clade or
#if it's not possible to have them all in one bipartition
def identify_ortho_issue(Folder,ingroup,logfile):

    array = []
    HASH = {}
    array = os.listdir(Folder)
    ingroup_array = []
    Extras.get_time("Dissecting Trees: ",logfile)
    
    outw = open("OutgroupOrthologyAnalysis.csv","w")
    outw2 = open("IngroupOrthologyAnalysis.csv","w")
    #get the gene conversion
    con_file = open("gene_conversion.csv","r")
    con_file.readline()
    for line in con_file:
        line = line.strip("\r\n")
        HASH[line.split(",")[1]] = line.split(",")[0]
    
    
    clades = []
    genome_seqs = []
    bipartition = []
    outw.write("Gene,WithGenome,OrthologyInvestigate,Genomes\n")
    outw2.write("Gene,WithGenome,OrthologyInvestigate,ReasonForSuggestion\n")
    for x in array:
        if x[6:14] == "bestTree":
            tree = open(Folder+x,"r")
            name_array,t = build(tree.readline())
            clades = post_order(t,name_array,clades)
            #how many sequences come from genomes
            genome_seqs = get_seqs_from_genomes(name_array)
            
            #get the smallest bipartition that contains all the genome
            #seqs in it
            bipartition,ortho_error = get_smallest_outgroup(clades,genome_seqs)
            gene = x.split(".")[1]
            outw.write(HASH[gene] + "," + gene + "," + str(ortho_error)+"\n")
            advice = False
            if ingroup != "":
                advice,ortho_error = get_smallest_ingroup_given(clades,genome_seqs,ingroup)
           
            outw2.write(HASH[gene] + "," + gene + "," + str(ortho_error) + "," + advice +"\n")
            
    outw.close()
    outw2.close()
#compare two arrays of trees in bipart form
def compare_array_biparts(array1,array2):

    total_biparts = len(array1)
    unique_to_tree1 = []
    
    for x in array1:
        l = x.split("|")[0]
        r = x.split("|")[1]
        l_array1 = l[:-1].split(",")
        r_array1 = r[1:].split(",")
        l_array1.sort()
        r_array1.sort()
        count = 0
        #run through the other array
        for y in array2:
            l = y.split("|")[0]
            l_array2 = l[:-1].split(",")
            l_array2.sort()
            if l_array1 == l_array2:
                count += 1
            elif r_array1 == l_array2:
                count += 1
        if count == 0:
            unique_to_tree1.append(x)
    return unique_to_tree1
            

#gets the clades that differ between the trees in question
def all_by_all_clades(array,logfile):
    
    logfile.write("Making a file called UniqueBipartitions.txt for species comp\n")
    outw = open("UniqueBipartitions.txt","w")
    edges = []
    Array = []
    for x in range(0,(len(array)-1)):
        for y in range((x+1),len(array)):
            temp_array = []
            temp_array.append("tree"+str(x)+"_tree"+str(y))
            unique1 = []
            unique2 = []
            unique1 = compare_array_biparts(array[x],array[y])
            temp_array.append(unique1[0])
            unique2 = compare_array_biparts(array[y],array[x])
            temp_array.append(unique2[0])
            outw.write("Comparison tree" + str(x) + "_tree" + str(y) + "\n")
            for z in unique1:
                outw.write("\tTree" + str(x+1) + ": " + str(z.replace(","," ") + "\n"))
            for z in unique2:
                outw.write("\tTree" + str(y+1) + ": " + str(z.replace(","," ") + "\n"))
            outw.write("\n")
            Array.append(temp_array)
    outw.close()
    return Array
            
#Get all the bipartitions of trees and find the ones unique to the trees
def get_dif_biparts(species_trees,logfile):
    
    count = 0
    count = Extras.get_lines(species_trees)
    file = open(species_trees,"r")
    aa = []
    Array = []
    
    for line in file:
        clades = []
        name_array = []
        
        line = line.strip("\r\n")
        name_array,t = build(line)
        clades = post_order(t,name_array,clades)

        aa.append(clades)
    Extras.get_time("Comparing all species trees to all",logfile)
    Array = all_by_all_clades(aa,logfile)
    file.close()
    return Array

#remove elements missing from one list and another
def rm_missing(species,name_array):
    
    return_array = []
    for x in species:
        if x in name_array:
            return_array.append(x)
    return return_array
        

#checks to see if out of a list of bipartitions
#the bipart in question exists
def is_bipart_in_list(bipart,list,name_array):

    Match = "False"
    l_species = bipart.split("|")[0][:-1].split(",")
    r_species = bipart.split("|")[1][1:].split(",")
    
    #account for missing data
    r_species_cor = rm_missing(r_species,name_array)
    l_species_cor = rm_missing(l_species,name_array)
    r_species_cor.sort()
    l_species_cor.sort()
    
    for x in list:
        l_gene = x.split("|")[0][:-1].split(",")
        l_gene.sort()
        #identify an identical bipartition
        if l_gene == r_species_cor or l_gene == l_species_cor:
            Match = "True"
    return Match

    



def get_trees_matching_species(sp_tree_bip,HASH_GWLL,tree_folder,logfile):
    
    #this will have the match between the GWLL file and the tree names in the folder
    numb = []
    
    logfile.write("Checking which gene trees match the species tree\n")
    logfile.write("Working to identify gene name to GWLL\n")
    numb = Extras.match_gene_name_to_GWLL(HASH_GWLL,tree_folder,logfile)
    
    HASH = {}
    array = os.listdir(tree_folder)
    for x in array:
        
        if x.split(".") == 0:
            name = x
        elif len(numb) == 1:
            name = x.split(".")[numb[0]]
        else:
            name = x.split(".")[numb[0]:numb[1]]
        if tree_folder[-1] == "/":
            file = open(tree_folder+x,"r")
        else:
            file = open(tree_folder+"/"+x,"r")
        name_array,t = build(file.readline())
        file.close()
        clades = []
        clades = post_order(t,name_array,clades)
        
        #check to see if a bipart in question exists in the gene tree
        t_count = 0
        for y in sp_tree_bip:
        
            match_gwll_tree1 = "True"
            match_gwll_tree2 = "False"
            match_tree1 = is_bipart_in_list(y[1],clades,name_array)
            match_tree2 = is_bipart_in_list(y[2],clades,name_array)
            if match_tree1 == "True" and match_tree2 == "True":
                match_tree1 = "Uninformative"
                match_tree2 = "Uninformative"
            if float(HASH_GWLL[name][t_count]) < 0.:
                match_gwll_tree1 = "False"
                match_gwll_tree2 = "True"

            if y[0] in HASH:
                HASH[y[0]].append(name + "," + match_gwll_tree1 + "," + match_gwll_tree2 + "," + match_tree1 + "," + match_tree2 + "," + HASH_GWLL[name][t_count] + "," + y[1].replace(","," ") + "," + y[2].replace(","," "))
            else:
                HASH[y[0]] = []
                HASH[y[0]].append(name + "," + match_gwll_tree1 + "," + match_gwll_tree2 + "," + match_tree1 + "," + match_tree2 + "," + HASH_GWLL[name][t_count] + "," + y[1].replace(","," ") + "," + y[2].replace(","," "))
            t_count += 1
    
    for x in HASH:
        tree1 = x.split("_")[0]
        tree2 = x.split("_")[1]
        outw = open("Comp_"+tree1+"_"+tree2+".csv","w")
        outw.write("GeneName,GWLLsupport_"+ tree1+",GWLLsupport_"+tree2+",MLsupport_"+tree1+",MLsupport_"+tree2+",GWLLvalue,bipartition_"+tree1+",bipartition_"+tree2+"\n")
        for y in HASH[x]:
            outw.write(y+"\n")
        outw.close()



