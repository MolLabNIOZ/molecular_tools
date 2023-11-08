# -*- coding: utf-8 -*-
"""
Automated shortening of sequences in FASTA format

VERSION: November2023

Dowload the asv sequences in a .txt or .seq file. Make sure the file is in 
FASTA format!
The name of the file should be the NIOZ number of the sequence lane.
This protocol uses the relative paths from the GitHub folder.
Add the relative path to this file behind 'pathway'.
If you want to use a file or save the file in a different folder, change 
the relative paths to a different folder.

edit:
231108: adds a repetition of 20 A's behind every sequence to hopefully make it 
easier to allign the sequences. Also added a function that allows you to 
disable the function that adds the poly A end after each sequence. It also adds
the domain name of the species.

@author: rdebeer
"""
#import the package named panda
import pandas

#you can customize the output of the protocol. 
extention = True #makes it possible to dissable the addition of the A-string, if False, it doesn't add the A-string. If true, it adds the A-string
tax = True #adds the domain of the species, if False, it adds unknown to the sequence

#checks if you want to add the domain name to the FASTA output
if tax:
#makes a variable named file for the file with the ASV table
    file = "molecular_tools/FASTA shortening/asvTable_noSingletons.txt"

#changes the .txt file to a csv file
    ASV_table = pandas.read_csv(file,delimiter="\t",header=1)

#creates a dictonary for the taxonomy per ASV number
    taxonomies = {} 

#adds the taxonomy to the dictonary per ASV number
    for i in ASV_table.index:
        ASV = ASV_table["#OTU ID"][i]
        taxonomy = (str((ASV_table["taxonomy"][i]))).split(";")[0]
        taxonomies[ASV] = taxonomy

#creates a variable for the file you want to use. You can use a .txt or .seq as starting point. the pathway to your file behind pathway
pathway = "molecular_tools/FASTA shortening/NIOZ354.seq"
#creates a variable for the NIOZ number using the name of the file in the pathway above
NIOZ = pathway[-11:-4]
#opens the file using the variable you created above and creates it under the variable 'data'
data = open(pathway)

#creates the output variable, this will be later used to export to the new .txt file
output1 = ""
AA = "" #a variable to add a 'n' amount of A's behind every sequence.


if extention: #executes the next lines if extention = True
    for i in range(0,20): #sets the amount of A's in the string
        AA = AA + 'A' #loop to create the string
    
#creates a loop that checks every line in the file you opened
for line in data:

#checks if the line stats with '>', if it does, it checks if the ASV number found on that line is the same as the ASV number in the dictonary.
#if this is not the case, it adds unknown as taxonomy
    if line.startswith('>'):
        try:    
            ASV = line[1:].strip() #.strip deletes the enter at the end of the line
            domain = taxonomies[ASV] 
        except:
            domain = "unknown" #adds unknown as domain to the sequence if the domain is not in the dictonary
        #generates the output you add for every line in the file to the variable you created in line 55    
        output1 = output1 + line.strip() + "_" + domain + "_" + NIOZ + "\n"
#if the line doesn't start with a >, it takes the first 50 characters on the line and adds this to the output. "\n" adds a new line to the output
    else:
        output1 = output1 +line[0:50] + AA + "\n"
    
#creates a new .seq file in the right folder. NIOZ is the variable you have created above.
#if you want to create a .txt file, change .seq to .txt
#'w' is a function for writing to the file. 'as f' creates a variable you use in the next line
#if you want to save the file somewhere else, change the first argument to another folder, leave the rest as it currently is.
with open('molecular_tools/FASTA shortening/verkorte data/' + NIOZ + "_shortened.seq", 'w') as f:
#writes the output variable you created in line 15 to the file you created
    f.write(output1)

#closes the .txt file so you can work with it outside Spyder
data.close()
f.close()
print("The shortening of your sequence data has been completed")