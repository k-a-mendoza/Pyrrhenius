import glob

import numpy as np
import pandas as pd
import os
import pickle
import matplotlib.pyplot as plt
file = 'Database/publication_database.csv'
from pyrrhenious import database
from pyrrhenious import utils as pyrutils

def get_fig_data(ec_model_row, image_directory):
    source_dir = os.sep.join([image_directory,str(ec_model_row['publication_id'].iloc[0]),'*.png'])
    pngs = glob.glob(source_dir)
    extents=[]
    aspects=[]
    xticks=[]
    yticks =[]
    for png in pngs:
        with open(png[:-4]+'_extent.pkl','rb') as f:
            data = pickle.load(f)
        extents.append(data['extent'])
        aspects.append((np.max(data['xticks'])-np.min(data['xticks']))/((np.max(data['yticks'])-np.min(data['yticks']))))
        xticks.append(data['xticks'])
        yticks.append(data['yticks'])
    return pngs, extents, aspects, xticks, yticks

def prep_fig(png, extent,xticks,yticks,aspect):
    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot()
    ax.imshow(plt.imread(png),extent=extent,aspect=aspect)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.grid(True)
    return fig, ax
