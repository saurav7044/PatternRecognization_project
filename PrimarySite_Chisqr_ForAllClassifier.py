#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
from sklearn.feature_selection import VarianceThreshold
from sklearn.naive_bayes import GaussianNB
from imblearn.over_sampling import RandomOverSampler, SMOTE, ADASYN
from collections import Counter
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis, LinearDiscriminantAnalysis


Clinical = pd.read_excel('C:/prad_tcga_clinical_data.xlsx')
GeneData = pd.read_excel('C:/prad_tcga_genes.xlsx',header=None)
GeneData = GeneData[~np.all(GeneData == 0, axis=1)]


# In[153]:


GeneTr = GeneData.T
Header = GeneTr.iloc[0]
GeneTr = GeneTr[1:]
GeneTr.columns = Header
CombinedData = pd.merge(GeneTr, Clinical, how="inner", left_on='ID', right_on='PATIENT_ID')
X = CombinedData.ix[:, 1:60484]
y = CombinedData['PRIMARY_SITE']


# In[154]:


#Feature Selection
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
clf = SelectKBest(score_func=chi2, k=100)     # DEFINE THE NUMBER OF FEATURES IN K
clf = clf.fit(X, y)

model = SelectFromModel(clf, prefit=True)


# In[155]:


X_new = clf.transform(X)


# In[156]:


#Over sampling
#ros = RandomOverSampler(random_state=0)
#X_new_1, y_new_1 = ros.fit_resample(X_new, y)
X_new_1, y_new_1 = SMOTE().fit_resample(X_new, y)
print(sorted(Counter(y_new_1).items()))
print("Shape of oversampled data " + str(X_new_1.shape))


# In[157]:


#Classification
from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis, LinearDiscriminantAnalysis
a = []
b = []
#INITIALIZE 25 ARRAYS TO STORE CONFUSE MATRIX CELL VALUES

#################      4    #############
def w1_w1(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Central Zone'):
            count+=1
    return count
    
def w1_w2(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Overlapping / Multiple Zones'):
            count+=1
    return count
    
def w1_w3(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Peripheral Zone'):
            count+=1
    return count
    
def w1_w4(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Transition Zone'):
            count+=1
    return count
    
def w1_w5(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == '[Not Available]'):
            count+=1
    return count


    
#################     2       #############
def w2_w1(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Central Zone'):
            count+=1
    return count
    
def w2_w2(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Overlapping / Multiple Zones'):
            count+=1
    return count
    
def w2_w3(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Peripheral Zone'):
            count+=1
    return count
    
def w2_w4(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Transition Zone'):
            count+=1
    return count
    
def w2_w5(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == '[Not Available]'):
            count+=1
    return count
        
#################     3      #############
def w3_w1(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Central Zone'):
            count+=1
    return count
    
def w3_w2(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Overlapping / Multiple Zones'):
            count+=1
    return count
    
def w3_w3(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Peripheral Zone'):
            count+=1
    return count
    
def w3_w4(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Transition Zone'):
            count+=1
    return count
    
def w3_w5(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == '[Not Available]'):
            count+=1
    return count


    
#################     4    #############
def w4_w1(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Central Zone'):
            count+=1
    return count
    
def w4_w2(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Overlapping / Multiple Zones'):
            count+=1
    return count
    
def w4_w3(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Peripheral Zone'):
            count+=1
    return count
    
def w4_w4(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Transition Zone'):
            count+=1
    return count
    
def w4_w5(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == '[Not Available]'):
            count+=1
    return count

#################     5      #############
def w5_w1(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Central Zone'):
            count+=1
    return count
    
def w5_w2(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Overlapping / Multiple Zones'):
            count+=1
    return count
    
def w5_w3(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Peripheral Zone'):
            count+=1
    return count
    
def w5_w4(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == 'Transition Zone'):
            count+=1
    return count
    
def w5_w5(y_true, y_pred,n,a):
    count=0
    for i in range(n):
        if (y_true[i] == a) and (y_pred[i] == '[Not Available]'):
            count+=1
    return count



W1W1 =[]
W1W2 =[]
W1W3 =[]
W1W4 =[]
W1W5 =[]

W2W1 =[]
W2W2 =[]
W2W3 =[]
W2W4 =[]
W2W5 =[]

W3W1 =[]
W3W2 =[]
W3W3 =[]
W3W4 =[]
W3W5 =[]

W4W1 =[]
W4W2 =[]
W4W3 =[]
W4W4 =[]
W4W5 =[]

W5W1 =[]
W5W2 =[]
W5W3 =[]
W5W4 =[]
W5W5 =[]
FinalConfuseMatrix = np.empty([5,5])
kf = KFold(n_splits=10,random_state=2,shuffle = True)
for train, test in kf.split(X_new_1):  # SPLIT TRAIN AND TEST DATA 
   # print("TRAIN:", train, "TEST:", test)
    X1, X2 = X_new_1[train], X_new_1[test]
    Y1, Y2 = y_new_1[train], y_new_1[test]
    y_true = Y2
    n= len(Y2)
    
    #clf = GaussianNB()
    clf = svm.SVC(kernel='rbf', C = 1.0, gamma='auto')
    #clf = svm.SVC(kernel='linear', C = 1.0)
    #clf  = KNeighborsClassifier(n_neighbors=3)
    #clf = QuadraticDiscriminantAnalysis()
    #clf= LinearDiscriminantAnalysis()
    #clf=RandomForestClassifier(n_estimators=50)
    y_pred = (clf.fit(X1,Y1).predict(X2))   
    
    w1w1 = w1_w1(y_true, y_pred,n,'Central Zone')
    w1w2 = w1_w2(y_true, y_pred,n,'Central Zone')
    w1w3 = w1_w3(y_true, y_pred,n,'Central Zone')
    w1w4 = w1_w4(y_true, y_pred,n,'Central Zone')
    w1w5 = w1_w5(y_true, y_pred,n,'Central Zone')
    
    w2w1 = w2_w1 (y_true, y_pred,n,'Overlapping / Multiple Zones')
    w2w2 = w2_w2 (y_true, y_pred,n,'Overlapping / Multiple Zones')
    w2w3 = w2_w3 (y_true, y_pred,n,'Overlapping / Multiple Zones')
    w2w4 = w2_w4 (y_true, y_pred,n,'Overlapping / Multiple Zones')
    w2w5 = w2_w5 (y_true, y_pred,n,'Overlapping / Multiple Zones')
    
    w3w1 = w3_w1 (y_true, y_pred,n,'Peripheral Zone')
    w3w2 = w3_w2 (y_true, y_pred,n,'Peripheral Zone')
    w3w3 = w3_w3 (y_true, y_pred,n,'Peripheral Zone')
    w3w4 = w3_w4 (y_true, y_pred,n,'Peripheral Zone')
    w3w5 = w3_w5 (y_true, y_pred,n,'Peripheral Zone')

    w4w1 = w4_w1 (y_true, y_pred,n,'Transition Zone')
    w4w2 = w4_w2 (y_true, y_pred,n,'Transition Zone')
    w4w3 = w4_w3 (y_true, y_pred,n,'Transition Zone')
    w4w4 = w4_w4 (y_true, y_pred,n,'Transition Zone')
    w4w5 = w4_w5 (y_true, y_pred,n,'Transition Zone')
    
    w5w1 = w5_w1 (y_true, y_pred,n,'[Not Available]')
    w5w2 = w5_w2 (y_true, y_pred,n,'[Not Available]')
    w5w3 = w5_w3 (y_true, y_pred,n,'[Not Available]')
    w5w4 = w5_w4 (y_true, y_pred,n,'[Not Available]')
    w5w5 = w5_w5 (y_true, y_pred,n,'[Not Available]')
    
    # APPEND THE VALUES FOR EACH ITERATION 
    W1W1.append(w1w1)
    W1W2.append(w1w2)
    W1W3.append(w1w3)
    W1W4.append(w1w4)
    W1W5.append(w1w5)

    W2W1.append(w2w1)
    W2W2.append(w2w2)
    W2W3.append(w2w3)
    W2W4.append(w2w4)
    W2W5.append(w2w5)

    W3W1.append(w3w1)
    W3W2.append(w3w2)
    W3W3.append(w3w3)
    W3W4.append(w3w4)
    W3W5.append(w3w5)

    W4W1.append(w4w1)
    W4W2.append(w4w2)
    W4W3.append(w4w3)
    W4W4.append(w4w4)
    W4W5.append(w4w5)

    W5W1.append(w5w1)
    W5W2.append(w5w2)
    W5W3.append(w5w3)
    W5W4.append(w5w4)
    W5W5.append(w5w5)
    
    print(str(accuracy_score(y_true, y_pred)*100) + "%")

# STORE THE FINAL CELL VALUES OF CONFUSE MATRIX
FinalConfuseMatrix[0][0] = np.sum(W1W1)
FinalConfuseMatrix[0][1] = np.sum(W1W2)
FinalConfuseMatrix[0][2] = np.sum(W1W3)                                    
FinalConfuseMatrix[0][3] = np.sum(W1W4)
FinalConfuseMatrix[0][4] = np.sum(W1W5)

FinalConfuseMatrix[1][0] = np.sum(W2W1)
FinalConfuseMatrix[1][1] = np.sum(W2W2)
FinalConfuseMatrix[1][2] = np.sum(W2W3)                                    
FinalConfuseMatrix[1][3] = np.sum(W2W4)
FinalConfuseMatrix[1][4] = np.sum(W2W5)
                  
FinalConfuseMatrix[2][0] = np.sum(W3W1)
FinalConfuseMatrix[2][1] = np.sum(W3W2)
FinalConfuseMatrix[2][2] = np.sum(W3W3)                                    
FinalConfuseMatrix[2][3] = np.sum(W3W4)
FinalConfuseMatrix[2][4] = np.sum(W3W5)

FinalConfuseMatrix[3][0] = np.sum(W4W1)
FinalConfuseMatrix[3][1] = np.sum(W4W2)
FinalConfuseMatrix[3][2] = np.sum(W4W3)                                    
FinalConfuseMatrix[3][3] = np.sum(W4W4)
FinalConfuseMatrix[3][4] = np.sum(W4W5)

FinalConfuseMatrix[4][0] = np.sum(W5W1)
FinalConfuseMatrix[4][1] = np.sum(W5W2)
FinalConfuseMatrix[4][2] = np.sum(W5W3)                                    
FinalConfuseMatrix[4][3] = np.sum(W5W4)
FinalConfuseMatrix[4][4] = np.sum(W5W5)
                  
                  
print("\nFinal Confuse Matrix of SVM\n",FinalConfuseMatrix)
#ConfuseMatrix = confusion_matrix(a, b)
#print(ConfuseMatrix)

#train the model
from sklearn.model_selection import cross_val_score

scores = []
scores = cross_val_score(clf,X1,Y1,cv=10)


# In[158]:


print(scores.mean())
print(str(scores.mean()*100) + "%")


# In[ ]:





# In[ ]:




