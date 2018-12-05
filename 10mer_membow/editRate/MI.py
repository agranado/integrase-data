
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import math
import numpy.matlib as matplot
import pandas as pa
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
            P_xy[x,y]=sum(np.logical_and(list1==alphabet[x], list2==alphabet[y]))/cell


    #Matrix method method from https://github.com/swainlab/mi-by-decoding/blob/master/matlab/MutualInformationCode/Info.m
    P_xy  = P_xy  + np.min(P_xy[P_xy>0]) *math.pow(10,-8) #for avoiding -inf values
    L = np.log2(P_xy)
    #entropy
    H = np.multiply(L,P_xy)

    Py = np.sum(P_xy,axis=0)
    Px = np.sum(P_xy,axis=1)

    Px_mat = matplot.repmat(Px,len(alphabet),1).transpose()
    Py_mat = matplot.repmat(Py,len(alphabet),1)
    #now the matrices have the right dimensions

    #Shannon's formula based on entropy: MI(X,Y) = H(X) - H(X|Y)
    I = sum(sum(H)) -  sum(sum(np.multiply(P_xy,  np.log2(np.multiply(Px_mat,Py_mat)))))
    return(I)




#MI NOTES:
#diagonal elements: The mutual information of X with itself is just its self information H(X)  = -P(X)*log(P(X)),
#https://math.stackexchange.com/questions/160621/what-is-the-mutual-information-ixx



filename = 'allBarcodes.csv'
df = pa.read_csv('allBarcodes.csv')
print(df)

#for big R just add a new character here
alphabet =['u','r','x','w']
colnames = list(df.columns.values)

MI_matrix = np.zeros((len(colnames),len(colnames)))

for i in range(0,len(colnames)):
    for j in range(0,len(colnames)):
        MI_matrix[i,j] = MI_xy( df[colnames[i]], df[colnames[j]],alphabet  )

#rename the data frame
MI_df.columns = pd.DataFrame(MI_mamtrix)
MI_df.index = colnames
MI_df.columns = colnames



# # # # # # # PLOT

fig,ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right',size='5%',pad=0.05)
im = ax.imshow(MI_df)
ax.set_title('MI all sites')
ax.set_xticks(np.arange(0,10,1))
ax.set_yticks(np.arange(0,10,1))
ax.set_yticklabels(colnames)
ax.set_xticklabels(colnames)
fig.colorbar(im,cax=cax,orientation ='vertical')
plot.show()

fig.savefig('MI_10mer.png', dpi = 1000)
