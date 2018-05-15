#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 10:26:59 2018

@author: chuck
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 14:54:28 2018

@author: chuck
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
#from sklearn.datasets import load_digits
import seaborn as sns; sns.set()
from sklearn import preprocessing

read_counts = pd.read_csv('/home/chuck/Documents/RNAseq/coding_norm/coding_gene_sums_norm.txt', sep='\t')
read_information = pd.read_csv('/home/chuck/Documents/RNAseq/RNAseq_sample_info/ltee_replicate1vsAll.txt', sep='\t')

features = read_counts['features']

replicate = read_information['replicate']

read_counts = read_counts.drop(['features'], 1)

for read in read_counts:
    read_counts[read] = read_counts[read].astype(int)

libs = read_information['sample'].values

colors = read_information['replicate'].values
read_counts_T = read_counts.T
read_counts_T.columns = features
    
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression

model = LogisticRegression()
rfe = RFE(model, 50)
fit = rfe.fit(read_counts_T, replicate)
print("Num Features: " + str(fit.n_features_)) 
print("Selected Features: " + str(fit.support_))
print("Feature Ranking: " + str(fit.ranking_))

out = pd.DataFrame([fit.ranking_, fit.support_], columns=read_counts_T.columns,index=['ranking', 'support']).T

out.to_csv('/home/chuck/Documents/RNAseq/coding_norm/RNAseq_feature_selection.txt', sep="\t", index=False)








#out = pd.DataFrame(pca.components_, columns=read_counts_T.columns,index=['PC1', 'PC2', 'PC3', 'PC4']).T


#Cen3D = plt.figure()
#ax = Cen3D.add_subplot(111, projection = '3d')
#
#ax.scatter(projected[:, 0], projected[:, 1], projected[:, 3], cmap='brg', c=colors)
#ax.set_xlabel('PC1')
#ax.set_ylabel('PC2')
#ax.set_zlabel('PC3')
#
#ax.view_init(25, -50)

