# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 11:41:28 2026

@author: rdebeer
"""
# =============================================================================
# Variables to fill in
# =============================================================================
# Link the absolute pathway to the .csv file
filepath = 'C:/Users/rdebeer/OneDrive - NIOZ/Data/python tests/nanodrop/Alg met Elsa.csv'

# Give your graph a title
graph_title = "Nanodrop van de 5 algen extracten"
# =============================================================================
# Imports the necessary 
import matplotlib.pyplot as plt
import pandas as pd
# Creates a dataframe from the .csv file with your data
dataframe = pd.read_csv(filepath)

# Creates a dataframe with only the data necessary for the linegraph (absorbance data)
dataframe_linegraph = dataframe.drop(dataframe.columns[[0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]], axis=1)

#### The procedure to create the plot
# Get the absorbance data and convert them to ints instead of strings
x = dataframe_linegraph.columns[1:].astype(int)

# Creates the figure
plt.figure(figsize=(10, 6))

# Loop for all the rows in the dataframe. This will create the lines per sample
for placeholder, row in dataframe_linegraph.iterrows():
    
    # y-values are the absorbance at a certain wavelenght
    y = row[1:].values
    
    # Plots the line for the sample
    plt.plot(x, y, alpha=0.6, label=row["Sample Name"])


# Creates the labels and graph title
plt.xlabel("Wavelength nm")
plt.ylabel("10 nm Absorbance")
plt.title(graph_title)

# Adds a legend, makes sure that the plat is as you want
plt.legend()
plt.xlim(220, 350)
plt.margins(x=0)

# Creates a grid on the plot
plt.gca().set_axisbelow(True)
plt.grid(True, linestyle="--", alpha=0.4)

# Creates a dict for the lines for the contaminations and DNA/RNA
lines = {
    230: ("Salts + organic compounds", "blue"),
    260: ("DNA/RNA", "green"),
    280: ("Proteins", "red")
}

# for-loop to add the 3 lines
for x_line, (label, color) in lines.items():
    
    plt.axvline(
        x=x_line,
        color=color,
        linestyle="--",
        linewidth=1,
        alpha=0.8
    )
    
    plt.text(
        x_line,
        plt.ylim()[1] * 0.95,
        label,
        rotation=90,
        va="top",
        ha="right",
        color=color
    )

plt.tight_layout()
plt.savefig(graph_title + '.png', dpi=300, bbox_inches="tight")
plt.show()

#%%
# Creates a dataframe with only the data for the infographic (this is without the absorbance data)
dataframe_infographic = dataframe.iloc[:, 1:20]