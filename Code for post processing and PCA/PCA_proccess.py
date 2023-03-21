# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 12:38:59 2022

@author: oscar
"""

## This PCA code was written to proccess the data generated for Automated Structural Activity Screening of β-Diketonate Assemblies with High-Throughput Ion Mobility - Mass Spectrometry
##imports
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

#settiungs for plotting
plt.rc('legend',fontsize=15)
font={'size':15,}
plt.rc('font',family='Arial')
plt.rcParams['axes.linewidth'] = 2


#read the relevant csvs
csv_1 = pd.read_csv("C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/csvs/Cu_Day10.csv") 
csv_2 = pd.read_csv("C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/csvs/Ni_Day12.csv") 
csv_3 = pd.read_csv("C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/csvs/Ni_Day3.csv") 
csv_4 = pd.read_csv("C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/csvs/Fe_Day7.csv")
csv_5 = pd.read_csv("C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/csvs/Cu_Day0.csv")
csv_6 = pd.read_csv("C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/csvs/Ni_M6_.csv")
csv_7 = pd.read_csv("C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/csvs/Fe_M6.csv")
csv_8 = pd.read_csv("C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/csvs/Cu_M6.csv")
csv_9 = pd.read_csv("C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/csvs/Fe_Day0.csv")
csv_10 = pd.read_csv("C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/csvs/Zn_Day0.csv")
csv_11=pd.read_csv("C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/csvs/Zn_Day3.csv")
csv_12=pd.read_csv("C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/csvs/Zn_Day10.csv")

#make them into a list 
combocsv=[csv_1,csv_2,csv_3,csv_4,csv_5,csv_6, csv_7,csv_8, csv_9,csv_10,csv_11,csv_12]
for csv in combocsv:
    f=[csv['deltaA'].max()]
    
    mai=csv.loc[csv['deltaA']==f[0]]
    g=mai['cluster']
    csv['cluster']=csv['cluster'].replace([int(g)],[4])

    i=[csv['meanmass'].max()]
    mai=csv.loc[csv['meanmass']==i[0]]
    g=mai['cluster']
    csv['cluster']=csv['cluster'].replace([int(g)],[6])

# making dataframe
combodf=pd.concat(combocsv)
combodf = combodf.rename(columns={'Unnamed: 0': 'expno'})

#removing very low deltaA
combodf.loc[combodf['deltaA']<0.8,'deltaA']=1

deltaaname= '\u0394A , Structural Diversity'
normalisedmassname=" Mean Normalised $\it{m/z}$"

combodf["metal"]=combodf["metal"].replace(["Cu","Ni", 'Fe','Zn'],[0,1,2,3])


##triming dataframe for pca
pcadf = combodf[['sample','metal','intensity', 'deltaA', 'Ratio','Ligandno','meanmass','time','pH']].copy()
pcadf = pcadf.reset_index(drop=True)


#measured variables for PCA
features = ['intensity', 'deltaA','meanmass']

# Separating out the features
x = pcadf.loc[:, features].values
# Separating out the target
y = pcadf.loc[:,['Ligandno']].values
# Standardizing the features
x = StandardScaler().fit_transform(x)

## PCA setup
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])

pca2 = PCA(n_components=2)
principalComponents2 = pca2.fit_transform(x)
principalDf2 = pd.DataFrame(data = principalComponents2
             , columns = ['principal component 1', 'principal component 2'])

finalDf = pd.concat([principalDf, pcadf[['Ligandno']]], axis = 1)
finalDf2 = pd.concat([principalDf2, pcadf[['metal']]], axis = 1)

#plot setup

plt.subplots(1,2,figsize=(20,8))
plt.subplot(121)
# rep = fig.add_subplot(1,1,1) 
plt.xlabel('PC1 (47.37 %)', fontsize = 20)
plt.ylabel('PC2 (33.95 %)', fontsize = 20)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)



# rep.set_zlabel('Principal Component 3', fontsize = 15)
#rep.set_title('2 component PCA', fontsize = 20)
targets = [ 2, 4,3,1,6,5]
legends=["L2", "L4", "L3","L1", "L5", "L6"]
colors = [ 'g', 'b','y','k','m','r']
plt.ylim(-2,6)

for target, color in zip(targets,colors):
    indicesToKeep = finalDf['Ligandno'] == target
    plt.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
                , finalDf.loc[indicesToKeep, 'principal component 2']
            
                , c = color
                , s = 50)
plt.legend(legends)

#fig2
plt.subplot(122)

plt.xlabel('PC1 (47.37 %)', fontsize = 20)
plt.ylabel('PC2 (33.95 %)', fontsize = 20)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
targets2 = [ 0,1,2,3]
legends2=["Cu","Ni", 'Fe','Zn']
colors2 = [ 'g', 'b','y','k']
plt.ylim(-2,6)
                

for target, color in zip(targets2,colors):
    indicesToKeep = finalDf2['metal'] == target
    plt.scatter(finalDf2.loc[indicesToKeep, 'principal component 1']
                , finalDf2.loc[indicesToKeep, 'principal component 2']
            
                , c = color
                , s = 50)
plt.legend(legends2)



#rep.grid()
print(pca.explained_variance_ratio_)
