import sys
import os
import Extras
import FastaTools

'''
Running the likelihood tests
'''
#gets the change in GWLL/SSLL
def get_delta_LL(array,name_array,outfile_name,logfile):

    Extras.get_time("Creating output for " + outfile_name.split(".")[0],logfile)
    outw = open(outfile_name,"w")
    
    tree = "GeneName,"
    #Create header
    for x in range(0,(len(array)-1)):
        for y in range((x+1),len(array)):
            tree += "tree" + str(x) + "_tree" + str(y) + ","
    outw.write(tree[:-1] + "\n")
    traverse = len(array[0])
    
    for x in range(0,traverse):
        cat = str(name_array[x]) + ","
        for y in range(0,(len(array)-1)):
            for z in range((y+1),len(array)):
                delta = 0.
                delta = float(array[y][x]) - float(array[z][x])
                cat += str(delta) + ","
        outw.write(cat[:-1] + "\n")
    Extras.get_time("Results written to " + outfile_name,logfile)
    outw.close()

#print out the data and summarize the results
def summarize_LL_arrays(array,name_array,outfile_name,logfile):

    Extras.get_time("Creating output for " + outfile_name.split(".")[0],logfile)
    tree = "GeneName,"
    outw = open(outfile_name,"w")
    
    #Create the header
    for x in range(0,(len(array)-1)):
        tree += "tree" + str(x) + ","
    tree += "tree" + str(x+1)
    outw.write(tree + "\n")

    #do an all by all comp of the GWLLs
    traverse = len(array[0])
    cat = ""
    
    for x in range(0,traverse):
        cat = name_array[x] + ","
        for y in range(0,len(array)):
            cat += str(array[y][x]) + ","
        outw.write(cat[:-1] + "\n")
    
    Extras.get_time("Results written to " + outfile_name,logfile)
    outw.close()

    

#Get the GWLL and SSLL from the site file
#and the parts
def get_GWLLand_SSLL(parts,site_file,logfile):

    part_array = []
    temp_array = []
    LL_array = []
    avgSSLL_array = []
    GWLL_array = []
    name_array = []

    
    first_line = True
    print "Getting the GWLL and Average SSLL info"
    Extras.get_time("Getting GWLL and Average SSLL",logfile)
    name_array,part_array = FastaTools.get_parts_array(parts)
    
    file = open(site_file,"r")
    for line in file:
        line = line.strip("\r\n")
        #ignore the first line of the site info
        if first_line == True:
            first_line = False
        else:
        
            LL_array = line.split()
            LL_array = [float(i) for i in LL_array[1:]]
            y = 0
            temp_GWLL_array = []
            temp_avgSSLL_array = []
            
            for x in part_array:
                
                z = x
                x += y
                temp_array = LL_array[y:x]
                GWLL = sum(temp_array)
                avgSSLL = (float(GWLL) / float(z))
                temp_GWLL_array.append(GWLL)
                temp_avgSSLL_array.append(avgSSLL)

            GWLL_array.append(temp_GWLL_array)
            avgSSLL_array.append(temp_avgSSLL_array)
    
    summarize_LL_arrays(GWLL_array,name_array,"GWLL.csv",logfile)
    summarize_LL_arrays(avgSSLL_array,name_array,"AvgSSLL.csv",logfile)
    
    get_delta_LL(GWLL_array,name_array,"deltaGWLL.csv",logfile)
    get_delta_LL(avgSSLL_array,name_array,"deltaAvgSSLL.csv",logfile)
    
    file.close()
    

#runs the two topology test
def RunRaxML_TwoTop(sp_tree,supermatrix,parts,raxml,logfile):

    Extras.get_time("Running two top",logfile)
    
    #get some basics on the matrix
    genes = 0
    taxa = 0
    trees = 0
    dna = True
    
    genes = Extras.get_lines(parts)
    taxa = Extras.get_lines(supermatrix,True)
    trees = Extras.get_lines(sp_tree)
    
    #check over the parts file and die if something is off
    p = open(parts,"r")
    dna = FastaTools.get_part_format(p.readline(),logfile)
    p.close()
    
    print "===================================================================="
    logfile.write("==================================================================\n")
    print "Analyzing " + str(genes) + " with " + str(taxa) + " taxa and " + str(trees) + " trees"
    logfile.write("Analyzing " + str(genes) + " with " + str(taxa) + " taxa and " + str(trees) + " trees\n")
    if dna == True:
        cmd = raxml + " -f G -T 2 -s " + supermatrix + " -q " + parts + " -m GTRGAMMA " + " -z " + sp_tree + " -n Topologies_SSLL"
    else:
        cmd = raxml + " -f G -T 2 -s " + supermatrix + " -q " + parts + " -m PROTGAMMAWAG " + " -z " + sp_tree + " -n Topologies_SSLL"
    print "Running your command: " + cmd
    logfile.write("Your raxml command " + cmd + "\n")
    #os.system(cmd)
    logfile.write("Command has run your results are in: RAxML_perSiteLLs.Topologies_SSLL\n")
    
    get_GWLLand_SSLL(parts,"RAxML_perSiteLLs.Topologies_SSLL",logfile)
    
    
def RunRaxml(Folder,raxml,dna,logfile):

    Extras.get_time("Starting Tree Inference",logfile)
    if dna == True:
        model = "GTRGAMMA"
    else:
        model = "PROTGAMMAWAG"
    #raxmlHPC-AVX -s TempOrthoFolder/gene1.fa.aln -p 12345 -m GTRGAMMA -n gene1
    array = []
    array = os.listdir(Folder)
    for x in array:
        array2 = x.split(".")
        if array2[-1] == "aln":
            cmd = raxml + "-T 2 -s " + Folder + x + " -p 12345 -m " + model + " -n " + array2[0]
            logfile.write("RaxML command: " + cmd + "\n")
            os.system(cmd)
            cmd = "mv " + "RAxML*" + array2[0] + " " + Folder
            os.system(cmd)
            

    
