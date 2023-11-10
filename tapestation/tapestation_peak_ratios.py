"""
This protocol will check the ratio between 16S and 18S bands based on
tapestation quantification after 515F/926R library prep

After doing a tapestation run, first make sure all lower markers are assigned 
correctly. If not, change the lower marker.
Then, in the Tapestation Analysis Software (4.1.1) press the [region] button 
(3rd button on the ribbon). Then right mouse click right on the 
electropherogram and press [Edit Region Settings]. Then press [Click here to 
add a new row]. Set the From and To [bs]s and choose a color. Then check the 
[Assay File] check box in order to apply the region to all samples in the file. 
Then press [Apply]. Press [Click here to add a new row] again to add the second
region.

Suggestion: 16S (300 - 520) / 18S (520 - 850)

Then press [File] --> [Export Data]. Uncheck all check boxes, except for
[Region Table] and  [Compact Region Table]. Then give your file a name and 
press [Export]. 
"""

# Import needed packages=======================================================
import pandas as pd
# =============================================================================

# Variables to set ============================================================
# How many regions do you want to analyze:
regions = 1 # 1 or 2
region_comment_1 = 'Euk'
if regions == 2:
    region_comment_2 = '18S'

# Where is the compact region table located:
filepath = ("//zeus.nioz.nl/mmb/molecular_ecology/mollab_team/Projects/2023/MMB/Linda/231011-NIOZ369_515F951R_Quant - 2023-10-11 - 13-49-06-D1000_compactRegionTable.csv")

# Data analysis ===============================================================
data = pd.read_csv(filepath, encoding='unicode-escape')

data_region_1 = (data[data['Region Comment'].
                      str.contains(region_comment_1)].
                     reset_index()
                     )
if regions == 2:
    data_region_2 = (data[data['Region Comment'].
                          str.contains(region_comment_2)].
                         reset_index()
                         )

df = pd.DataFrame({'Sample Description': data['Sample Description'].unique()})

df[region_comment_1 + '_ng/µL'] = data_region_1['Conc. [ng/µl]']
df[region_comment_1 +'band'] = ''

if regions == 2:
    df[region_comment_2 + '_ng/µL'] = data_region_2['Conc. [ng/µl]']
    df[region_comment_2 +'band'] = ''

for sample in df.index:
    conc_region_1 = df[region_comment_1 + '_ng/µL'][sample]
    
    
    if conc_region_1 < 0.2:
        band_region_1 = "no or very faint band"
    elif 0.2 < conc_region_1 < 1.5:
        band_region_1 = "faint band"
    elif 1.5 < conc_region_1 < 2:
        band_region_1 = "band"
    else:
        band_region_1 = "clear band"
    
    df.at[sample,region_comment_1 +'band'] = band_region_1
    
    if regions == 2:
        conc_region_2 = df[region_comment_2 + '_ng/µL'][sample]
        if conc_region_2 < 0.2:
            band_region_2 = "no or very faint band"
        elif 0.2 < conc_region_2 < 1.5:
            band_region_2 = "faint band"
        elif 1.5 < conc_region_2 < 2:
            band_region_2 = "band"
        else:
            band_region_2 = "clear band"
            
        
        if 0.8 < conc_region_2 / conc_region_1 < 1.5:
            result = region_comment_1 + ' = ' + region_comment_2
        
        else: 
            if conc_region_1 > conc_region_2:
                result = region_comment_1 + ' > ' + region_comment_2
            else:
                result = region_comment_1 + ' < ' + region_comment_2

        df.at[sample,region_comment_1 + region_comment_2 +' ratio'] = result
        df.at[sample,region_comment_2 +'band'] = band_region_2
    
## Save results
# Save file as region_comment_1 + region_comment_2
if regions == 2:
    filename = region_comment_1 + '_' + region_comment_2
else:
    filename = region_comment_1

df.to_csv(filename + '_bands.csv', sep="\t", index=False)  
