#make heatmap using correlation matrix

from string import ascii_letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from MI_barcodes import MI_xy
import math



##
filenames = ['output_all.csv','output_no_r.csv','output_no_x.csv']
for f in range(0,len(filenames)):
    df = pd.read_csv(filenames[f])
    #print(df)

    # alphabet =['u','r','x','w']
    alphabet = pd.unique( df.values.ravel() )
    colnames = list(df.columns.values)

    MI_matrix = np.zeros((len(colnames),len(colnames)))
    ### this is the same for all codes
    for i in range(0,len(colnames)):
        for j in range(0,len(colnames)):
            MI_matrix[i,j] = MI_xy( df[colnames[i]], df[colnames[j]],alphabet  )
    #here we have the matrix
    MI_df= pd.DataFrame(MI_matrix)
    MI_df.index = colnames
    MI_df.columns = colnames
    ###
    #plot from: https://seaborn.pydata.org/examples/many_pairwise_correlations.html
    sns.set(style="white")


    # Generate a mask for the upper triangle
    mask = np.zeros_like(MI_df, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_title('MI ' + filenames[f])
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    # cmap = sns.palplot(sns.cubehelix_palette(8),as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    max_val = MI_matrix[np.where(~np.eye(MI_df.shape[0],dtype=bool))].max()
    max_val = math.ceil(max_val *100)/100.0

    sns.heatmap(MI_df, mask=mask, cmap=cmap, vmax=max_val, center=0,
                square=True, linewidths=.5, cbar_kws={"shrink": .5})

    fig.savefig(filenames[f]+'_MI_OFFdiag.png', dpi = 1000)
