from datetime import datetime
import sys
import os

#Basic time thingy
def get_time(position, outfile):
	now = datetime.now()
	dt_string = now.strftime("%B %d, %Y %H:%M:%S")
	if outfile:
		print position + " " + dt_string
		outfile.write(position + " " + dt_string + "\n")
	else:
		print position + " " + dt_string

#Find out how many lines a file has
def get_lines(file,fasta=False):
    
    count = 0
    file1 = open(file,"r")
    for line in file1:
        if fasta == False:
            count += 1
        else:
            if line[0] == '>':
                count += 1
    file1.close()
    return count

#Return a hash of the delta GWLL
def get_deltaGWLL(file_1):

    HASH = {}
    file = open(file_1,"r")
    array = []
    
    for line in file:
        line = line.strip("\r\n")
        array = line.split(",")
        if len(array[0].split(".")) != 0:
            HASH[array[0].split(".")[0]] = array[1:]
        else:
            HASH[array[0]] = array[1:]
        
    return HASH
    
def match_gene_name_to_GWLL(HASH_GWLL,Folder,logfile):

    get_time("Working to identify how names match",logfile)
    test_array = []
    array = []
    numb_array = []
    array = os.listdir(Folder)
    
    if array[0].split(".") == 0:
        if array[0] in HASH_GWLL:
            print "This appears to be the match between GWLL file and tree names: " + array[0]
            logfile.write("This appears to be the match between GWLL file and tree names: " + array[0] + "\n")
            numb_array.append(-1)
            return numb_array
    
    test_array = array[0].split(".")
    
    for x in range(0,(len(test_array)-1)):
        if test_array[x] in HASH_GWLL:
            numb_array.append(x)
            print "This appears to be the match between GWLL file and tree names: " + test_array[x]
            logfile.write("This appears to be the match between GWLL file and tree names: " + test_array[x] + "\n")
            return numb_array
            
    print "Looking for a more complex match"
    for x in range(0,(len(test_array)-1)):
        for y in range((x+1),len(test_array)):
            if ".".join(test_array[x:y]) in HASH_GWLL:
                print "This appears to be the match between GWLL file and tree names: " + ".".join(test_array[x:y])
                logfile.write("This appears to be the match between GWLL file and tree names: " + ".".join(test_array[x:y]) + "\n")
                numb_array.append(x)
                numb_array.append(y)
                return numb_array
    print "No match could be found, the program examines what is delimited by a . so the gene tree files need to have the name within them"
    logfile.write("No match could be found, the program examines what is delimited by a . so the gene tree files need to have the name within them\n")
    sys.exit(0)

def check_if_happened(Item,message):
    
    isFile = os.path.exists(Item)
    if isFile == True:
        print message
    return isFile

def check_if_happened_in_dir(Folder,pattern,pos,message):

    array = os.listdir(Folder)
    isFile = False
    for x in array:
        
        if pos == "end":
            test = len(x) - len(pattern)
            if x[test:] == pattern:
                isFile = True
        if pos == "beg":
            if x[:len(pattern)] == pattern:
                isFile = True
    if isFile == True:
        print message
    return isFile
    
