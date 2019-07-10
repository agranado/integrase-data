# uncompyle6 version 3.2.6
# Python bytecode 3.6 (3379)
# Decompiled from: Python 3.5.5 | packaged by conda-forge | (default, Jul 23 2018, 23:45:11)
# [GCC 4.2.1 Compatible Apple LLVM 6.1.0 (clang-602.0.53)]
# Embedded file name: /Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer_membow/editRate/MI_barcodes.py
# Compiled at: 2018-12-05 11:34:24
# Size of source mod 2**32: 3520 bytes
import matplotlib.pyplot as plt, scipy.stats as stats, numpy as np, math, numpy.matlib as matplot, seaborn as sns, pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap

#this function applies a Normalization by a constant = log2(3)
#the result units, therefore are trits (max 1 for integrase data)
def MI_xy(list1, list2, alphabet):
    P_xy = np.zeros((len(alphabet), len(alphabet)))
    MI = 0
    N = len(list1)
    for x in range(0, len(alphabet)):
        for y in range(0, len(alphabet)):
            P_xy[(x, y)] = sum(np.logical_and(list1 == alphabet[x], list2 == alphabet[y])) / N

    P_xy = P_xy + np.min(P_xy[P_xy > 0]) * math.pow(10, -8)
    L = np.log2(P_xy)/np.log2(3)
    H = np.multiply(P_xy, L)
    Py = np.sum(P_xy, axis=0)
    Px = np.sum(P_xy, axis=1)
    Px_mat = matplot.repmat(Px, len(alphabet), 1).transpose()
    Py_mat = matplot.repmat(Py, len(alphabet), 1)
    I = sum(sum(H)) - sum(sum(np.multiply(P_xy, np.log2(np.multiply(Px_mat, Py_mat))/np.log2(3))))
    return I


#filenames = [
# 'output_all.csv', 'output_no_r.csv', 'output_no_x.csv']
filenames = ['allBarcodes.csv','output_all.csv']
filenames = ['output_all.csv']
filenames = ['allBarcodes.csv']
n_bootstrap = 100
n_sample = 500
b =0
for f in range(0, 1):
    df = pd.read_csv(filenames[f])
    #let's downsample the data frame with n_sample cells
    #so every time we do the analysis we do it with n_sample random cells and we do it n_bootstrap times
    df = df.sample(n_sample,axis = 0)
    alphabet = pd.unique(df.values.ravel())
    colnames = list(df.columns.values)
    MI_matrix = np.zeros((len(colnames), len(colnames), n_bootstrap))
    for b in range(0,n_bootstrap):
        for i in range(0, len(colnames)):
            for j in range(0, len(colnames)):
                MI_matrix[(i, j, b)] = MI_xy(df[colnames[i]], df[colnames[j]], alphabet)

MI_matrix_mean = np.mean(MI_matrix,axis = 2)
MI_df = pd.DataFrame(MI_matrix_mean)
MI_df.index = colnames
MI_df.columns = colnames
fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
#sns.cubhelix(10) return the defaul (pink) gradient
#start and rot determine other colors
#https://seaborn.pydata.org/tutorial/color_palettes.html
# ListedColormap converts a palette to a cmap to imshow can recognize it, otherwise you get error
#im = ax.imshow(MI_df,vmin = 0, vmax = 1, aspect='auto',cmap = ListedColormap(sns.cubehelix_palette(100))) #with scaling to max min values
# blue/green palette
im = ax.imshow(MI_df,vmin = 0, vmax = 1, aspect='auto',cmap = ListedColormap(sns.cubehelix_palette(100,start=.5, rot=-.75))) #with scaling to max min values

#im = ax.imshow(MI_df) #auto scaling (to see the actual range)
#sns.set_palette("Reds")
#im.set_cmap(sns.cubehelix_palette(10))


ax.set_title('MI ' + filenames[f])
ax.set_xticks(np.arange(0, 10, 1))
ax.set_yticks(np.arange(0, 10, 1))
ax.set_yticklabels(colnames)
ax.set_xticklabels(colnames)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.show()
fig.savefig(filenames[f] + '_MI_10mer.eps', dpi=300, format='eps')
# okay decompiling MI_barcodes.cpython-36.pyc
