#%% import module for pdf files
import pdfx

#%% Read pdf and create dictionary and list
# Read pdf file
pdf = pdfx.PDFx(
    "//ZEUS/mmb/molecular_ecology/mollab_team/Sequencing/ngs_sequencing/"
    "project_administration/NIOZ325/NIOZ325_RawData_report_HN00169497.pdf"
    )

# Get urls from pdf file as a dictionary
links_dict = pdf.get_references_as_dict()

# Convert dictionary to list
links_list = list(links_dict.values())

# Check what is printed to the file
# Because for some reason we get a list inside a list we need to use
# links_list[0] for the list.
for element in links_list[0]:
    print(element)
    
#%% Transfer list to txt file
#Create new txt file
txt_file = open(
    "download_links.txt", 
    "w"
    )
#Write list to txt with enter after each line
for element in links_list[0]:
    txt_file.write(element + '\n')
#close txt
txt_file.close()
