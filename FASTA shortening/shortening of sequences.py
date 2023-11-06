# -*- coding: utf-8 -*-
"""
Automated shortening of sequences in FASTA format

VERSION: November2023

Dowload the asv sequences in a .txt or .seq file. Make sure the file is in FASTA format!
The name of the file should be the NIOZ number of the sequence lane.
This protocol uses the relative paths from the GitHub folder.
Add the relative path to this file behind 'pathway'.
If you want to use a file or save the file in a different folder, change the relative paths to a different folder.

edit:


@author: rdebeer
"""
#creates a variable for the file you want to use. You can use a .txt or .seq as starting point. the pathway to your file behind pathway
pathway = "molecular_tools/FASTA shortening/NIOZ354.seq"
#creates a variable for the NIOZ number using the name of the file in the pathway above
NIOZ = pathway[-11:-4]
#opens the file using the variable you created above and creates it under the variable 'data'
data = open(pathway)

#creates the output variable, this will be later used to export to the new .txt file
output1 = ""

#creates a loop that checks every line in the file you opened
for line in data:
#if the line starts with '>', it adds the whole line to the output and adds _ and the NIOZ number with the variable NIOZ.
#line.strip() strips the layout of the line because this will generate the right output. "\n" adds a new line to the output
    if line.startswith('>'):
        output1 = output1 + line.strip() + "_" + NIOZ + "\n"
#if the line doesn't start with a >, it takes the first 50 characters on the line and adds this to the output. "\n" adds a new line to the output
    else:
        output1 = output1 +line[0:50] + "\n"
    
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