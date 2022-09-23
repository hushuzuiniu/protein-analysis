#!/bin/usr/python3

import sys, subprocess, shutil, os
import string, re
import collections
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import OrderedDict

my_ws = os.getcwd()     #to get the path of working space

################################# FUNCTIONS ################################
#before printing error trps I print ERROR message in order to get attention of the user
def error_msg(): 			
	print('\n')
	for i in range(20):        #twenty times print **  before word, want no new line in the end, because need to print word
		print("**",end="")
	print("ERROR", end="")      #then print word ERROR without newline
	for i in range(20):         #twenty times print ** after word
		print("**", end="")  
	print("\n")     #newline in the end          
	return (0)     #function just returns 0


def progress_msg():
        print('\n')
        for i in range(20):        #twenty times print ~~  before word, want no new line in the end, because need to print word
                print("~~",end="")
        print("PROGRESS", end="")      #then print word Progress without newline
        for i in range(20):         #twenty times print ~~ after word
                print("~~", end="")
        print("\n")     #newline in the end
        return (0)     #function just returns 0

#in order to check if the input is valid name, I check that word in easearch and count the found data. This check separately done for both inputs		
def check_valid(name):     #input is any name
	#when you do esearch it print out ENTREZ_DIRECT, where you can find found data count, I try to save that value and check
	count = subprocess.check_output('esearch -db protein -query "{0}" | xtract -pattern ENTREZ_DIRECT -element Count'.format(name),shell=True) 
	if int(count) == 0:      #if there is no data found thenfunctionss return 0 further it gaves error message to the user
		return 0
	return 1         #if there is at least one data found regarding to the user input we continue program
 
#now use both valid input names to check if there is data with that organism and with that protein names
def check_both_input(name,protein):   
	#saving as a value the Count data from esearch commamd output 
        count = subprocess.check_output('esearch -db protein -query "{0}[prot]" -organism "{1}" | xtract -pattern ENTREZ_DIRECT -element Count'.format(protein, name),shell=True)
        if int(count) == 0:       #if count 0 for input searched it means do data in the server with those inputs
                return 0	#to check function output it must return 0 when count 0
        return 1		#otherwise it return 1 for positive results

#this function for checking syntax errors and special characters then runs to esearch
def check_input(name):       #as an input we need one name
	valid = re.compile("[A-Za-z0-9- ]+")       #valid names can contain Upper or lower letters and numbers and space and dashes
	if name.isdigit() == True or valid.fullmatch(name) == None or re.search(r'^[\s\n\t]', name):      #error if there is only numbers or other special charaters or starts with sapce or tabs
		error_msg()     #calling error messaging  printing function
		print("The name "+name+" is invalid. Special characters, Only numbers, Whitespaces in the begining strongly restricted!\n")   #Explaining error for wrong input
		return 0     	#after finding error stops funciton and returns 0
	if check_valid(name) != 1:       #after checking for syntaxes, we call function which checks for esearch for valid data for that name, if it's not 1 then means count is 0
		error_msg()    #printing error msg          
		print("In database there is no imformation about "+name+". The name might have mispelling. Please try again!\n")  #Explaining the user the error   
		return 0      #after finding error stops funciton and returns 0  
	return 1  	#if all checking is correct then function returns 1

#if in the begining the inputs did not pass the checkings program ask for valid name
def get_input():     
	name = str(input("Please enter the taxonomic group name (No need for quotation!):\n"))    #asking to type the organism name
	while check_input(name) != 1: #until we dont get right input we call for function that checks syntax and in esearch combined
		name = str(input("Enter the valid taxonomic group name (No need for quotation):\n")) #if every input is wrong ask again for right input 
	protein = str(input("Please specify protein family name:\n")) #asking for protein name
	while check_input(protein) != 1:	  #until we dont get right input we call for function that checks syntax and in esearch combined	
		protein = str(input("Enter the valid protein name\n"))    #if every input is wrong ask again for right input
	return (name, protein)      #in the it returns the valid name for organism and protein

#before downloading all datas we run esearch with valid names and count all sequences and species after ask user if this what s/he wants
def begin(name, protein):   #need organism and protein name
	while check_input(name) != 1: #until we dont get right input we call for function that checks syntax and in esearch combined
		name = str(input("Enter the valid taxonomic group name (No need for quotation):\n")) #if every input is wrong ask again for right input
	while check_input(protein) != 1:          #until we dont get right input we call for function that checks syntax and in esearch combined
		protein = str(input("Enter the valid protein name\n"))    #if every input is wrong ask again for right input
	while check_both_input(name, protein) != 1:   #calling the function which ckeck for both names in esearch, until there is data for both input
		error_msg()		#calling function to print Error
		#Explaining the error message
		print("No data found with a taxonomic group name \""+name+"\" and protein family name \""+protein+"\" in database. Try one more time with the valid names!\n")
		name, protein = get_input()     #again going back, where we ask for inputs and check, then send here
	rename = name.replace(" ","_")     #if there is spaces in name I replace them to underscore in order to use it later for kaming file names
	reprotein = protein.replace(" ","_")  
	os.mkdir('{0}/{1}_{2}'.format(my_ws,rename, reprotein))   #By names with spaces replaced to undersocre, I make a folder for this run 
	os.chdir('{0}/{1}_{2}'.format(my_ws,rename, reprotein))  #Move to that folder, in order make it working directory
	#From all found Data about proteins sequenced for that taxonomic groups, first downloaded the organism list and save into file
	subprocess.call('esearch -db protein -query "{0}[prot]" -organism "{1}" | efetch -db protein -format docsum | xtract -pattern DocumentSummary -element Organism > {2}_{3}_organism_list.txt'.format(protein, name, rename,reprotein), shell=True)
	my_list = open('{0}_{1}_organism_list.txt'.format(rename, reprotein)).readlines()  #open organism list file and read by line saving it into array 
	count = len(my_list)    #length of read lines is organism numbers found
	species = len(set(my_list))  #founding unique in that list is for finding the number of species
	print('The search with your input resulted in '+str(count)+' sequences from '+str(species)+' species.\n')  #Telling the user about how many species are represented in the dataset chosen by user 
	if count <= 2:    #this program does not allow to use the dataset with more than 10,000 sequences
		error_msg()    #calling Error printer function
		print("Program will not for for less than 2 sequences!\n")
		reply = "BREAK"     #in this error case the function will return as BREAK
		return (reply, name, protein)    #it does not continue the program returns values
	if count >= 10000:    #this program does not allow to use the dataset with more than 10,000 sequences
		error_msg()    #calling Error printer function
		print("Allowable starting sequence set shouldn't be bigger than 10,000 sequences.\n")  
		reply = "BREAK"     #in this error case the function will return as BREAK
		return (reply, name, protein)    #it does not continue the program returns values
	reply = input('Do you want to continue the process with this dataset? (Otherwise you can stop it and run again with other inputs) yes/no\n').upper()  #otherwise program asks if user is satisfied with found dataset
	return (reply, name, protein)   #the user reply and organism and protein names transferred as output

############################### START #######################################################
#starting message
for i in range(40):    
	print('@@',end="")
print('\n\n')
print("PROGRAM IS STARTING TO RUN\n")

arguments = len(sys.argv) - 1    #taking out the program name
if (arguments != 2):      #to check if user inputted two names only
        error_msg()      #calling error printer function
        print("You have to enter two arguments to run the program!Quote the Names which have spaces!\n")
        sys.exit("Number of arguments are not valid! Please read the manual carefully! Program is exiting!\n")       #exiting the program until user enters two names for search

name  = sys.argv[1]     #fisrt argument after the program name is taxonomic group name
protein = sys.argv[2]	#second argument after the program name is protein name
	
reply, name, protein = begin(name, protein)     #calling my function that checks and if valid downloads organism list and returns the user reply and valid organism and protein name
while not (re.search('[YES]|[Y]',reply)):      #if the user reply is other then Yes or Y, we are returning to the beginning 
	print("We are going back to the beginning!")
	rename = name.replace(" ","_")         #Since we dont return renamed names from previous functions, here we are doing it again  
	reprotein = protein.replace(" ","_")   #replacing spaces with underscores
	os.chdir('{0}'.format(my_ws))	       #going out of the folder which was created for this run
	shutil.rmtree('{0}/{1}_{2}'.format(my_ws,rename, reprotein))   #since user does not want to use this dataset we are deleting whole directory
	name, protein = get_input()            #getting new inputs
	reply, name, protein = begin(name,protein)  #calling function begin with new outputs, and inside while loop will check again for reply is yes gets out of the loop

#Okay now starting to Download fasta files
inside_dir = os.getcwd()	 #setting path of the folder with tax_group and protein name
result_dir = ("{0}/RESULTS".format(inside_dir))     #all significant results will be saved in new folder called RESULTS
os.mkdir(result_dir)
progress_msg()
print("All sequences are being downloaded...")      
rename = name.replace(" ","_")    #Since we dont return renamed names from previous functions, here we are doing it again
reprotein = protein.replace(" ","_")    #replacing spaces with underscores
#downloading fasta file with that orgainsm and protein name
subprocess.call('esearch -db protein -query "{0}[prot]" -organism "{1}" | efetch -db protein -format fasta > {2}_{3}.fasta'.format(protein, name, rename, reprotein), shell=True)    
#using clustalo to align all sequences
subprocess.call('clustalo -i {0}_{1}.fasta -o {0}_{1}_aln.fasta -v'.format(rename, reprotein), shell =True)
progress_msg()
print("Multiple Alignment for downloaded protein sequences is done\n")

#Asking user if want to see the alignment result in pretty format
show_al = str(input("Do you want to visualize alignment results in pretty format? (y/n)\n"))
if (re.search('[YES]|[Y]',show_al.upper())):    #for positive results we run prettyplot functin for aligned sequences
	os.chdir(result_dir)  #in order to save plot output changing to RESULT folder
	subprocess.call('prettyplot -residuesperline=70 -boxcol -consensus -ratio=0.59 -docolour -sequence {0}/{1}_{2}_aln.fasta -graph svg'.format(inside_dir,rename,reprotein),shell=True)        #in each line wanted to 70 residues and, respresent residues and backgrounds in colors, for consensus match used 0.59 as pluraity ratio
	progress_msg()
	print("Program is drawing the aligned sequences with pretty formating!\n")	
	subprocess.call('firefox {0}/prettyplot.svg&'.format(result_dir), shell=True)   #need to open output in firefox since plot is large
	os.chdir(inside_dir)

#For conservation plotting we are using the most similar 250 sequences from dataset
#Running cons function from 
# to create the consensus sequence from multiple alignment 
subprocess.call('cons -plurality 0.59 -sequence {0}_{1}_aln.fasta -outseq {0}_{1}_consensus.fasta -auto'.format(rename, reprotein), shell=True)   #used same prularity as for pretty plot in order to have same consensus sequence and it should ne a bit more than half of the total sequence weighting 0.59
progress_msg()
print("Program has created a consensus sequence from a multiple alignment\n\n")

#We can get most similar sequences by running the blastp
#before running we are creating a database from all sequences
subprocess.call('makeblastdb -in {0}_{1}.fasta -dbtype prot -out {0}_{1}_db'.format(rename, reprotein), shell=True)
progress_msg()
print("From all fasta sequences the program created database for blast alingment\n")
#After, we are running the blast program with database and created consensus sequence
subprocess.call('blastp -db {0}_{1}_db -query {0}_{1}_consensus.fasta -outfmt 7 -max_hsps 1 > {0}_{1}_similarity_seq_blast.out'.format(rename, reprotein), shell = True)   #outformat wanted with tabs and comments, per each query and subject sequences we want only one HSPs 
progress_msg()
print("Alighned all sequences in BLAST and their HSP score saved in new file\n")

#opening the output of the blast file
blast_file = open("{0}/{1}_{2}_similarity_seq_blast.out".format(inside_dir,rename, reprotein)).read().rstrip('\n')
access_hsp = {}    #settyng empty dictionary
data_lines = blast_file.split('\n')  #splitting file into lines 
for lines in data_lines:            #for each line
	if re.search('#',lines):     # if comment skip the line
		next
	else:
		data_tab = lines.split('\t')      #split the line columns by tab-separator
		acc = data_tab[1]                 #second column is accession numbers, we are saving into value
		hsp = data_tab[11]        #last column is HSP scores
		access_hsp[acc] = hsp     #for each accession numbers saving the score into dictionary
access_hsp_ord = {}         #will save the ordered 250 most similar accession numbers and score in different dictionary
#reversely, in descending order [until 250] sorting values of the dictionary and saving to new dict
access_hsp_ord=OrderedDict(sorted(access_hsp.items(), key=lambda value: value[1], reverse=True)[:250])
acc_list_file = open("acc_filt_header.txt", "w")   #opening new file to save the output of the sorted 250 sequence accession number
acc_list_file.write("\n".join(access_hsp_ord.keys())) #ordered keys written by new line to the new file
acc_list_file.close()  #closing file
#created file with 250 accession numbers used to pull matching sequeces from single long fasta file
#here using the pullseq program which located in Assignment2 file in BPSM directory
subprocess.call('/localdisk/data/BPSM/Assignment2/pullseq -i {0}_{1}_aln.fasta -n acc_filt_header.txt -v > similar_aln_seq_250.fasta'.format(rename, reprotein), shell=True) #to give extra details during run used flag -v

#insert window size for plotting
progress_msg()
print("The most similar (maximum 250) sequences were pulled from the fasta file. Program will plot the conservation level!\n")
wind_size = str(input("Please enter the WINDOW size for plotting the level of conservation: (any integer) \n"))
wind_valid = re.compile("[0-9]+")   #checking if user eneter only digits
if wind_valid.fullmatch(wind_size) == None:     #if not then set default value 
	wind_size = 10
	print("Invalid window size, we will use default size (10)!\n")     

#plotting the conservation level of most similar sequences using EMBOSS plotcon function
os.chdir(result_dir)
subprocess.call('plotcon -sequence {0}/similar_aln_seq_250.fasta -winsize {1} -graph svg'.format(inside_dir,wind_size), shell=True)   
subprocess.call('eog {0}/plotcon.svg&'.format(result_dir), shell = True)    #'eog' used for presenting output to the screen and '&' sign used so that program continues its funtion further, while action of presenting figure taking place in the behind of main program
os.chdir(inside_dir)
progress_msg()
print("Plotcon result went to the screen and saved in RESULTS folder!\n")

#for calling my files I need to put accession numbers as a name for each files
fasta_all = open("{0}_{1}.fasta".format(rename, reprotein)).read().rstrip()   #open single fasta file 
os.mkdir('{0}/FASTA_FILES'.format(inside_dir)) 			#making directory for each separate fasta files
os.chdir('{0}/FASTA_FILES'.format(inside_dir))			#moving to new directory place
pat_out = ('{0}/PATMOTIFS_OUT'.format(inside_dir))		#for PATMOTIFS program also creating new directory
os.mkdir(pat_out)
found_motif = open('{0}/FOUND_MOTIFS.txt'.format(result_dir),"w")  #creating file to save found MOTIF information	
found_motif.write("Accession number\tMotif name\n")		 #The header of the file will print names of columns
motif_array = []   	 #empty array
each_seq = fasta_all.split(">") 	#spliting all fasta files by char >
for seq in each_seq:			#for each sequence in fasta file
	for keys in access_hsp.keys():	 #for each accession number in dictionaries from blast output
		if re.match(keys, seq):	  #checking for match in order to tag each sequence with accession number
			if re.search('|', keys):  #strange symbols in accession numbers are replaced by space while creating file
				keys = keys.replace("|","_")
			sep_file = open("{0}.fasta".format(keys), "w")	#tagged accesion numbers will be used to name each file for sequence
			sep_file.write(">"+seq)				#after creating file we write seqences into file, while splittin char '>' disappeared, I am putting it back 
			sep_file.close()	#after writing closing the file
			#each written sequence file used to run PATMATMOTIFS from EMboss to scan the protein sequence for known motifs
			subprocess.call("patmatmotifs -sequence {0}.fasta -outfile {1}/{0}.patmatmotifs -full".format(keys, pat_out), shell=True) #used -full flag to run for whole sequence
			pat_file = open("{0}/{1}.patmatmotifs".format(pat_out, keys)).readlines() #after we run the patmotif we open output to check for motif
			for line in pat_file:		#through each lines
				if re.search('#',line):		#skipping comments
					next
				elif re.search('Motif', line):	     #searching for line where Motif word used
					line.rstrip()	             #cutting new line from the end of line, in order to not to get in my match	
					index = line.find("=")       #motif name is after the sign '=', that's why I want to get index of this sign 
					motif = line[index+2:]	     #getting name of the Motif after '= ' and space, so I add 2 
					motif_array.append(motif)    #saving found motif names into array		
					found_motif.write('{0}\t{1}'.format(keys,motif))  #samely, writing into the file together with accession number
found_motif.close()    #after writing to file need to close
os.chdir(inside_dir)  #after change of folders coming back to working space
progress_msg()
print("The output of all sequences' scanning for motifs from the PROSITE database are saved in the folder by name PATMOTIFS_OUT inside the  working directory\n")
print("In subset of protein sequences {1} known motifs were found associated and unique motif is {0}. In FOUND_MOTIFS.txt file you can find list of accession numbers with found motif names. Repeated name means that protein sequence has more than one times the same motif.\n".format(len(set(motif_array)),len(motif_array)))
print("Found MOTIF names are: ")
#found Motif information is presented to user
print("\n".join(set(motif_array)))

#asking if user want to analyze the secondary structure of the protein sequences

sec_confirm = str(input("Do you analyze the Secondary srtucture of the protein sequences? Percentege of the secondary structure will be presented (y/n)\n"))
if (re.search('[YES]|[Y]',sec_confirm.upper())): #if yes then run the program
	sec_str = ('{0}/Secondary_STRUCT'.format(inside_dir))  #for Secondary structure output creating new directory 
	os.mkdir(sec_str)
	H = {}      #empty dictionaries for each structure form Helix
	E = {}		#Sheet
	T = {}		#Turns
	C = {}		#Coils
	for keys in access_hsp.keys():	#for each fasta files
		if re.search('|', keys):
			keys = keys.replace("|","_")
		#running Emboss garnier function for each sequence to predict the secondary structure
		subprocess.call("garnier -sequence {0}/FASTA_FILES/{1}.fasta -outfile {2}/{1}.garnier".format(inside_dir,keys,sec_str), shell=True)
	second_files = os.listdir('{0}'.format(sec_str))  #all output files saved into one list 
	for item in second_files:		#for each file name in Secondary_STRUCT folder
		keys = item.replace(".garnier","")    #first get only accession numbers by cutting the end part of the name
		sec_file = open("{0}/{1}".format(sec_str,item)).read().rstrip()    #open each file to read
		#from file look for line where The percentage for each structure represented
		m = re.search(r"percent: H: (-?\d+.\d+) E: (-?\d+.\d+) T: (-?\d+.\d+) C: (-?\d+.\d+)",sec_file) 
		#each matched percentage value saved according to their Structure names in separate dictionaries by their accession number		  
		H[keys] = m.group(1)
		E[keys] = m.group(2)
		T[keys] = m.group(3)
		C[keys] = m.group(4)
	#making new text file to write the statistics of the each struct
	progress_msg()
	print("The Program predicted the secondary structures for all the fasta sequences!\n")
	table = open("{0}/Secondary_struct_percent_table.txt".format(result_dir), "w")  
	#writing the header line for file
	table.write("Accession numbers\tHelix %\tSheet %\tTurns %\tcoils %\n")
	#for each accession number writing each dictionary value by tab-separator 
	for key in H.keys():
       		table.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(key,H[key],E[key],T[key],C[key]))
	table.close()  
	print("The statistical summary about second structure for each sequence saved in RESULTS folder.\n")
	print("Soon plot will appear in your screen! To end the program please close the plot window!")

	#reading as a dataframe the table that contains percentages of for each sequences
	S = pd.read_csv("{0}/Secondary_struct_percent_table.txt".format(result_dir), sep = "\t") 
	S.plot(figsize=(12,12))    #plotting the datafram by setting the figure size
	plt.xlabel('Sequence numbers')    #name of x axis
	plt.ylabel('Percentage')	 #name for y axis
	plt.title('Secondary structural elements distribution in Protein sequence dataset')  #title of the graph
	plt.savefig("{0}/Secondary_struct_profile_plot.png".format(result_dir), transparent=True)  #saving the figure in the Result folder
	plt.show()   #putting image out


print('\n')
for i in range(20):        #twenty times print ##  before word, want no new line in the end, because need to print word
	print("##",end="")
print("END", end="")      #then print word Progress without newline
for i in range(20):         #twenty times print ## after word
	print("##", end="")
print("\n")     #newline in the end

print('This is the end of the PROGRAM! Thank you for using! Hope to see you again\n')
