# -*- coding: utf-8 -*-
"""

"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ABC_toolbox import ABC_plot, ABC_utils, cell_funcs

#%%

local_data_dir = r'C:\Users\jpv88\Documents\GenePicker9001_data'
freqs = pd.read_pickle(os.path.join(local_data_dir, "antIRN-PARN-MERFISH-freqs.pkl"))

#%%

class_freqs = cell_funcs.recalc_freqs(freqs, new_level='class')
subclass_freqs = cell_funcs.recalc_freqs(freqs, new_level='subclass')
supertype_freqs = cell_funcs.recalc_freqs(freqs, new_level='supertype')

# %%

filepath = os.path.join("MapMyCells",
                        "MapMyCellsInput_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1777211802468",
                        "MapMyCellsInput_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1777211802468.csv")

MapMyCellsResults = pd.read_csv(filepath, comment='#')

MapMyCellsSubclassFreqs = cell_funcs.calc_frac_per_type(MapMyCellsResults, level="subclass_name")


# %%

sub_ratios = {k: v for k, v in sorted(subclass_freqs.items(), key=lambda item: item[0])}
sub_color_dict = ABC_plot.fetch_colors_dict("subclass")
colors = [sub_color_dict[x] for x in sub_ratios.keys()]
    
x = np.array(list(sub_ratios.values()))

fig, ax = plt.subplots()
ax.pie(x, colors=colors, wedgeprops=dict(width=0.25),startangle=90)

# %%

sub_ratios = {k: v for k, v in sorted(MapMyCellsSubclassFreqs.items(), key=lambda item: item[0])}
sub_color_dict = ABC_plot.fetch_colors_dict("subclass")
colors = [sub_color_dict[x] for x in sub_ratios.keys()]
    
x = np.array(list(sub_ratios.values()))

fig, ax = plt.subplots()
ax.pie(x, colors=colors, wedgeprops=dict(width=0.25),startangle=90)

# %%

MapMyCellsSubclass = MapMyCellsResults["subclass_name"].values
real_subclasses = list(subclass_freqs.keys())

inIRNPARN = 0
outIRNPARN = 0

for subclass in MapMyCellsSubclass:
    if subclass in real_subclasses:
        inIRNPARN += 1
    else:
        outIRNPARN += 1
        
n_cells = len(MapMyCellsResults)

# %%

import matplotlib as mpl

labels = ["inside ant-IRN-PARN", "outside ant-IRN-PARN"]

bars = plt.bar(labels,[inIRNPARN, outIRNPARN])

for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, height, str(height), ha="center", va="bottom")
    
plt.ylabel("# cells")
plt.title("MapMyCells ABC subclasses results")

mpl.rcParams['image.composite_image'] = False
plt.rcParams['svg.fonttype'] = 'none'

#%%

# Sample data
categories = real_subclasses
group_a = list(subclass_freqs.values())
group_b = []

for real_subclass in real_subclasses:
    if real_subclass in MapMyCellsSubclassFreqs:
        group_b.append(MapMyCellsSubclassFreqs[real_subclass])
    else:
        group_b.append(0)
        
group_a = np.array(group_a) * 100
group_b = np.array(group_b) * 100
        
x = np.arange(len(categories))
width = 0.35  # width of each bar
 
fig, ax = plt.subplots(figsize=(8, 5))
 
bars_a = ax.bar(x - width/2, group_a, width, label='ABC Atlas Cell Type Proportions', color='steelblue')
bars_b = ax.bar(x + width/2, group_b, width, label='MapMyCells Cell Type Proportions', color='coral')

# Labels and formatting
ax.set_xlabel('ABC Atlas ant-IRN-PARN Subclasses')
ax.set_ylabel('Cell Type Proportion')
ax.set_title('ABC Atlas vs. MapMyCells')
ax.legend()

