
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import math
import numpy.matlib as matplot
import seaborn as sns

import pandas as pd

from mpl_toolkits.axes_grid1 import make_axes_locatable



def MI_xy(list1,list2,alphabet):

    P_xy = np.zeros((len(alphabet),len(alphabet)))
    MI = 0

    #for each pair of columns we need P(x,y) the join distribution which is a matrix
    #of probabilities
    #P(x,y) = P(X=x, Y=y) for all values, in this case x,u,r,w


    N = len(list1)
    for x in range(0, len(alphabet)):
        for y in range(0, len(alphabet)):
            P_xy[x,y]=sum(np.logical_and(list1==alphabet[x], list2==alphabet[y]))/N


    #Matrix method method from https://github.com/swainlab/mi-by-decoding/blob/master/matlab/MutualInformationCode/Info.m
    P_xy  = P_xy  + np.min(P_xy[P_xy>0]) *math.pow(10,-8) #for avoiding -inf values
    L = np.log2(P_xy)
    #entropy
    H = np.multiply(P_xy,L)

    Py = np.sum(P_xy,axis=0)
    Px = np.sum(P_xy,axis=1)

    Px_mat = matplot.repmat(Px,len(alphabet),1).transpose()
    Py_mat = matplot.repmat(Py,len(alphabet),1)
    #now the matrices have the right dimensions

    #Shannon's formula based on entropy
    I = sum(sum(H)) -  sum(sum(np.multiply(P_xy,  np.log2(np.multiply(Px_mat,Py_mat)))))
    return(I)




#MI NOTES:
#diagonal elements: The mutual information of X with itself is just its self information H(X)  = -P(X)*log(P(X)),
#https://math.stackexchange.com/questions/160621/what-is-the-mutual-information-ixx



filenames = ['output_all.csv','output_no_r.csv','output_no_x.csv']

#for f in range(0,len(filenames)):
for f in range(0,1):

    df = pd.read_csv(filenames[f])
    #print(df)

    # alphabet =['u','r','x','w']
    alphabet = pd.unique( df.values.ravel() )
    colnames = list(df.columns.values)

    MI_matrix = np.zeros((len(colnames),len(colnames)))

    for i in range(0,len(colnames)):
        for j in range(0,len(colnames)):
            MI_matrix[i,j] = MI_xy( df[colnames[i]], df[colnames[j]],alphabet  )


    MI_df= pd.DataFrame(MI_matrix)
    MI_df.index = colnames
    MI_df.columns = colnames


    # # # # # # #

    fig,ax = plt.subplots()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right',size='5%',pad=0.05)
    im = ax.imshow(MI_df)
    ax.set_title('MI ' + filenames[f])
    ax.set_xticks(np.arange(0,10,1))
    ax.set_yticks(np.arange(0,10,1))
    ax.set_yticklabels(colnames)
    ax.set_xticklabels(colnames)
    fig.colorbar(im,cax=cax,orientation ='vertical')
    plt.show()

    fig.savefig(filenames[f]+'_MI_10mer.png', dpi = 1000)
    # matplotlib.pyplot.close()

#######

#plot heatmap
#sns.heatmap(MI_df,annot=True)


# fig1 = plt.figure() # create a figure with the default size
#
#
# ax1 = fig1.add_subplot(1,2,1)
# im = ax1.imshow(MI_df, interpolation='none')
#
# ax1.set_title('MI all sites')
#
# ax1.set_xticks(np.arange(0,10,1))
# ax1.set_yticks(np.arange(0,10,1))
# ax1.set_yticklabels(colnames)
# ax1.set_xticklabels(colnames)

#ax2 = fig1.add_subplot(2,2,2)
#ax2.imshow(im2, interpolation='none')
#ax2.set_title('100 X 100')




#Number of effective states given the probabilities
#if it was 2 state distribution otherwise use log3, log4
# np.power(2,np.diagonal(MI_matrix))

#python Notes:
#Select all rows and some cols:
# df.loc[:,'1':'10']  //  '1':'10' is a slice // could be a list df.loc[:,['1','6','10']]

#for frequency we can also use groupby
#aa=df.groupby('2').size()
 #then access using aa['u'] or aa[alphabet[0]]
