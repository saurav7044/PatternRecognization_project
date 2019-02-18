#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.datasets import load_iris
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from imblearn.over_sampling import RandomOverSampler
from sklearn import decomposition


# In[2]:


Clinical = pd.read_excel('C:/prad_tcga_clinical_data.xlsx')
GeneData = pd.read_excel('C:/prad_tcga_genes.xlsx',header=None)
GeneData = GeneData[~np.all(GeneData == 0, axis=1)]


# In[14]:


GeneTr= GeneData.T
Header = GeneTr.iloc[0]
GeneTr = GeneTr[1:]
GeneTr.columns = Header
GeneTr.drop(["GLEASON_PATTERN_PRIMARY", "GLEASON_PATTERN_SECONDARY", "GLEASON_SCORE", "CLIN_T_STAGE", "PATH_T_STAGE", "Transcript"], axis=1, inplace=True)
CombinedData = pd.merge(GeneTr, Clinical, how="inner", left_on='ID', right_on='PATIENT_ID')


# In[15]:


X = CombinedData.ix[:,1:60484]
y = CombinedData['PRIMARY_SITE']
print (y)


# In[16]:


ros = RandomOverSampler(random_state=0)
X, y = ros.fit_resample(X, y)


# In[17]:


#fig = plt.figure(1, figsize=(8, 8))
fig = plt.figure(1, figsize=(8, 8))
plt.clf()
ax = Axes3D(fig, rect=[0, 0, .95, 1],azim=60, elev=15)
plt.cla()
pca = decomposition.PCA(n_components=3)
pca.fit(X)
X = pca.transform(X)
for i in range(len(y)):
    if(y[i] == 'Overlapping / Multiple Zones'):
        ax.scatter3D(X[i][1], X[i][2],X[i][0], c='red',cmap=plt.cm.nipy_spectral)
    elif(y[i]=='Peripheral Zone'):
        ax.scatter3D(X[i][1], X[i][2],X[i][0], c='blue',cmap=plt.cm.nipy_spectral)
    elif(y[i] == 'Central Zone'):
        ax.scatter3D(X[i][1], X[i][2],X[i][0], c='green',cmap=plt.cm.nipy_spectral)
    elif(y[i] == 'Transition Zone'):
        ax.scatter3D(X[i][1], X[i][2],X[i][0], c='Orange',cmap=plt.cm.nipy_spectral)
    else:
        ax.scatter3D(X[i][1], X[i][2],X[i][0], c='pink',cmap=plt.cm.nipy_spectral)


plt.show()


# In[ ]:




