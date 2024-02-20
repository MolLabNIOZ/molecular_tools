# -*- coding: utf-8 -*-
"""
A protocol for converting the excel file with the students, users and PI's to
an output with all the current lab-users of the molecular lab.

You can choose which people you want to include or exclude in your email by 
changing the variables.
@author: rdebeer
"""

# Importing panda module, you need this for opening a xlsx and converting it to a dataframe
import pandas as pd

# Creating the different output based on the people you want to receive the mail
users = True      # Adds all the PhD'ers and PostDocs that are working at the MolLab
PI = True         # Adds all the PI's that are supervising any people working at the MolLab
students = True   # Adds the students that are currently working at the MolLab

# The variable for the file you want to use
xlsx = '//zeus/mmb/molecular_ecology/mollab_team/Administratie/Userinfo/Gebruikers, studenten, PI s en projecten.xlsx'

# Creates an empty list where all the generates mailadresses will be added to.
receivers = []

# If you want to add the users, it reads the xlsx file at the sheet named "Users" in colum B (the names) and it creates a dataframe
# Then it changes the dataframe to a list after the "Who"
if users is True:
    user = pd.read_excel(xlsx, 'Users', usecols="B")
    list_user = user["Who"].tolist()
    # For all the items in the list, it will replace all the spaces to a . and then adds @nioz.nl to the end and adds it to the empty list "receivers".
    for i in list_user:
        receivers.append(i.replace(" ", ".") + "@nioz.nl")
        
# If you want to add the users, it reads the xlsx file at the sheet named "PI" in colum B (the names) and it creates a dataframe
# Then it changes the dataframe to a list after the "Who"  
if PI is True:
    prof = pd.read_excel(xlsx, 'PIs', usecols="B")
    list_prof = prof["Who"].tolist()
    # For all the items in the list, it will replace all the spaces to a . and then adds @nioz.nl to the end and adds it to the empty list "receivers".
    for i in list_prof:
        receivers.append(i.replace(" ", ".") + "@nioz.nl")
        
# If you want to add the users, it reads the xlsx file at the sheet named "Students" in colum B (the names) and it creates a dataframe
# Then it changes the dataframe to a list after the "Who"
if students is True:
    student = pd.read_excel(xlsx, 'Students', usecols="B")
    list_student = student["Who"].tolist()
    # For all the items in the list, it will replace all the spaces to a . and then adds @nioz.nl to the end and adds it to the empty list "receivers".
    for i in list_student:
        receivers.append(i.replace(" ", ".") + "@nioz.nl")

# Prints the output receivers and changes the output of the list (changes ' ', to a ;)
print("; ".join(receivers))
print(len(receivers))