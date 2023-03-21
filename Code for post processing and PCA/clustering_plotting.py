# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 10:06:19 2021

@author: oscar
"""
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import AgglomerativeClustering

plt.rc('legend',fontsize=15)
font={'size':15,'family':'sans-serif','sans-serif':['Arial']}


plt.rc('font',**font)
plt.rcParams['axes.linewidth'] = 2

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

#
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

combodf=pd.concat(combocsv)
combodf = combodf.rename(columns={'Unnamed: 0': 'expno'})
# combodf['deltaA'] = combodf['deltaA'].apply(lambda x: 1 if x < 0.8)
combodf.loc[combodf['deltaA']<0.8,'deltaA']=1

mlcols=combodf[['meanmass','intensity','deltaA']]
newcl=mlcols.copy()

cluster2 = AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
cluster2.fit_predict(newcl)

combodf['newclus']=cluster2.labels_

deltaaname= '\u0394A , Structural Diversity'
normalisedmassname=" Mean Normalised $\it{m/z}$"

# combodf[]

#********************************************************************************#

x1 = combodf[ combodf['Ligandno'] == 1].deltaA.iloc[:]
y1 = combodf[ combodf['Ligandno'] == 1].meanmass.iloc[:]
z1 = combodf[ combodf['Ligandno'] == 1].intensity.iloc[:]
x2 = combodf[ combodf['Ligandno'] == 2].deltaA.iloc[:]
y2 = combodf[ combodf['Ligandno'] == 2].meanmass.iloc[:]
z2 = combodf[ combodf['Ligandno'] == 2].intensity.iloc[:]
x3 = combodf[ combodf['Ligandno'] == 3].deltaA.iloc[:]
y3 = combodf[ combodf['Ligandno'] == 3].meanmass.iloc[:]
z3 = combodf[ combodf['Ligandno'] == 3].intensity.iloc[:]
x4 = combodf[ combodf['Ligandno'] == 4].deltaA.iloc[:]
y4 = combodf[ combodf['Ligandno'] == 4].meanmass.iloc[:]
z4 = combodf[ combodf['Ligandno'] == 4].intensity.iloc[:]
x5 = combodf[ combodf['Ligandno'] == 5].deltaA.iloc[:]
y5 = combodf[ combodf['Ligandno'] == 5].meanmass.iloc[:]
z5 = combodf[ combodf['Ligandno'] == 5].intensity.iloc[:]
x6 = combodf[ combodf['Ligandno'] == 6].deltaA.iloc[:]
y6 = combodf[ combodf['Ligandno'] == 6].meanmass.iloc[:]
z6 = combodf[ combodf['Ligandno'] == 6].intensity.iloc[:]



shapesagain2=plt.figure()
axagain2=Axes3D(shapesagain2)
c=axagain2.scatter(x3,y3,z3, c='purple',label='L1')
a=axagain2.scatter(x1,y1,z1, c='y', label= 'L2', marker='D')
d=axagain2.scatter(x4,y4,z4, c='tomato',label='L3',marker='h')
b=axagain2.scatter(x2,y2,z2, c='k',label='L4',marker='s')
f=axagain2.scatter(x6,y6,z6, c='springgreen',label='L5', marker='P')
e=axagain2.scatter(x5,y5,z5, c='b',label='L6',marker="*")

axagain2.set_xlabel(deltaaname,fontname='Arial',fontsize=25, labelpad=12) #delatA
axagain2.set_ylabel(normalisedmassname,fontname='Arial',fontsize=25, labelpad=12) #maxmass
axagain2.set_zlabel('$Log_{10}$ Intensity',fontname='Arial',fontsize=25, labelpad=12) #intensity
axagain2.set_zlim3d(5,8)
axagain2.set_ylim3d(1,3.5)
axagain2.set_xlim3d(0.75,4.5)
axagain2.zaxis.set_tick_params(labelsize=15)
axagain2.yaxis.set_tick_params(labelsize=15)
axagain2.xaxis.set_tick_params(labelsize=15)

handles, labels = axagain2.get_legend_handles_labels()
# totallen=len(x1)+len(x2)+len(x3)+len(x4)+len(x5)+len(x6)
# patch = mpatches.Patch(color='k', label=(str(totallen)+' Spectra'))
# handles.append(patch)
axagain2.legend(handles=handles)

axagain2.view_init(elev=35, azim=145)
plt.show()

#********************************************************************************#

aaverage ,baverage,caverage= np.mean(x1),np.mean(y1),np.mean(z1)
daverage ,eaverage,faverage= np.mean(x2),np.mean(y2),np.mean(z2)
gaverage ,haverage,iaverage= np.mean(x3),np.mean(y3),np.mean(z3)
javerage ,kaverage,laverage= np.mean(x4),np.mean(y4),np.mean(z4)
maverage ,naverage,oaverage= np.mean(x5),np.mean(y5),np.mean(z5)
paverage ,qaverage,raverage= np.mean(x6),np.mean(y6),np.mean(z6)


shapesagain6=plt.figure()
axagain6=Axes3D(shapesagain6)
c=axagain6.scatter(gaverage ,haverage,iaverage,  c='purple',label='L1 Average',s=60)
a=axagain6.scatter(aaverage ,baverage,caverage, c='y', label= 'L2 Average', marker='D', s=60)
d=axagain6.scatter(javerage ,kaverage,laverage, c='tomato',label='L3 Average',marker='h',s=60)
b=axagain6.scatter(daverage ,eaverage,faverage, c='k',label='L4 Average',marker='s',s=60)
f=axagain6.scatter(paverage ,qaverage,raverage, c='springgreen',label='L5 Average', marker='P',s=60)
e=axagain6.scatter(maverage ,naverage,oaverage, c='b',label='L6 Average',marker="*",s=60)

axagain6.set_xlabel(deltaaname,fontsize=15) #delatA
axagain6.set_ylabel(normalisedmassname,fontsize=15) #maxmass
axagain6.set_zlabel('$Log_{10}$ Intensity',fontsize=15) #intensity

axagain6.set_zlim3d(5,8)
axagain6.set_ylim3d(1,4)
axagain6.set_xlim3d(0.75,4.5)

handles, labels = axagain6.get_legend_handles_labels()
totallen=len(x1)+len(x2)+len(x3)+len(x4)+len(x5)+len(x6)
patch = mpatches.Patch(color='k', label=(str(totallen)+' Spectra'))
handles.append(patch)
axagain6.legend(handles=handles)
axagain6.view_init(elev=35, azim=145)


plt.show()


#********************************************************************************#

combodf["metal"]=combodf["metal"].replace(["Cu","Ni", 'Fe','Zn'],[0,1,2,3])

x1 = combodf[ combodf['metal'] == 0].deltaA.iloc[:]
y1 = combodf[ combodf['metal'] == 0].meanmass.iloc[:]
z1 = combodf[ combodf['metal'] == 0].intensity.iloc[:]
x2 = combodf[ combodf['metal'] == 1].deltaA.iloc[:]
y2 = combodf[ combodf['metal'] == 1].meanmass.iloc[:]
z2 = combodf[ combodf['metal'] == 1].intensity.iloc[:]
x3 = combodf[ combodf['metal'] == 2].deltaA.iloc[:]
y3 = combodf[ combodf['metal'] == 2].meanmass.iloc[:]
z3 = combodf[ combodf['metal'] == 2].intensity.iloc[:]
x4 = combodf[ combodf['metal'] == 3].deltaA.iloc[:]
y4 = combodf[ combodf['metal'] == 3].meanmass.iloc[:]
z4 = combodf[ combodf['metal'] == 3].intensity.iloc[:]

cs=combodf['metal']
shapesagain3=plt.figure(figsize=(10,10))
axagain3=Axes3D(shapesagain3)
a=axagain3.scatter(x1,y1,z1, c='b', label= 'Cu', marker='D')
b=axagain3.scatter(x2,y2,z2, c='yellowgreen', label= 'Ni', marker='s')
c=axagain3.scatter(x3,y3,z3, c='lightcoral', label= 'Fe', marker='P')
d=axagain3.scatter(x4,y4,z4, c='black', label= 'Zn', marker='h')

axagain3.set_xlabel(deltaaname,fontname='Arial',fontsize=25, labelpad=12) #delatA
axagain3.set_ylabel(normalisedmassname,fontname='Arial',fontsize=25, labelpad=12) #maxmass
axagain3.set_zlabel('$Log_{10}$ Intensity',fontname='Arial',fontsize=25, labelpad=12) #intensity
axagain3.set_zlim3d(5,8)
axagain3.set_ylim3d(1,3.5)
axagain3.set_xlim3d(0.75,4.5)
axagain3.zaxis.set_tick_params(labelsize=15)
axagain3.yaxis.set_tick_params(labelsize=15)
axagain3.xaxis.set_tick_params(labelsize=15)

axagain3.legend(prop={'size': 18})
axagain3.view_init(elev=35, azim=145)
plt.show()


#********************************************************************************#
aaverage ,baverage,caverage= np.mean(x1),np.mean(y1),np.mean(z1)
daverage ,eaverage,faverage= np.mean(x2),np.mean(y2),np.mean(z2)
gaverage ,haverage,iaverage= np.mean(x3),np.mean(y3),np.mean(z3)
javerage ,kaverage,laverage= np.mean(x4),np.mean(y4),np.mean(z4)

shapesagain8=plt.figure(figsize=(10,10))
axagain8=Axes3D(shapesagain8)
a=axagain8.scatter(aaverage ,baverage,caverage, c='b', label= 'Cu Average', marker='D',s=60)
b=axagain8.scatter(daverage ,eaverage,faverage, c='yellowgreen', label= 'Ni Average', marker='s',s=60)
c=axagain8.scatter(gaverage ,haverage,iaverage, c='lightcoral', label= 'Fe Average', marker='P',s=60)
c=axagain8.scatter(javerage ,kaverage,laverage, c='black', label= 'Zn Average', marker='h',s=60)

axagain8.set_xlabel(deltaaname,fontname='Arial',fontsize=25, labelpad=12) #delatA
axagain8.set_ylabel(normalisedmassname,fontname='Arial',fontsize=25, labelpad=12) #maxmass
axagain8.set_zlabel('$Log_{10}$ Intensity',fontname='Arial',fontsize=25, labelpad=12) #intensity
axagain8.set_zlim3d(5,8)
axagain8.set_ylim3d(1,3.5)
axagain8.set_xlim3d(0.75,4.5)
axagain8.zaxis.set_tick_params(labelsize=15)
axagain8.yaxis.set_tick_params(labelsize=15)
axagain8.xaxis.set_tick_params(labelsize=15)

axagain8.legend( prop={'size': 18})
axagain8.view_init(elev=35, azim=145)

plt.show()

#********************************************************************************#
#pH
x1 = combodf[ combodf['pH']==3].deltaA.iloc[:]
y1 = combodf[ combodf['pH'] == 3].meanmass.iloc[:]
z1 = combodf[ combodf['pH'] == 3].intensity.iloc[:]
x2 = combodf[ combodf['pH']==7].deltaA.iloc[:]
y2 = combodf[ combodf['pH'] == 7].meanmass.iloc[:]
z2 = combodf[ combodf['pH'] == 7].intensity.iloc[:]
x3 = combodf[ combodf['pH']==10].deltaA.iloc[:]
y3 = combodf[ combodf['pH'] == 10].meanmass.iloc[:]
z3 = combodf[ combodf['pH'] == 10].intensity.iloc[:]

shapesagain4=plt.figure()
axagain4=Axes3D(shapesagain4)
lo=axagain4.scatter(x1,y1,z1, c='r', label= 'pH 3', marker='D')
mid=axagain4.scatter(x2,y2,z2, c='g',label='pH 7',marker='s')
hi=axagain4.scatter(x3,y3,z3, c='purple',label='pH 10')

axagain4.set_xlabel(deltaaname) #delatA
axagain4.set_ylabel(normalisedmassname) #maxmass
axagain4.set_zlabel('$Log_{10}$ Intensity') #intensity
axagain4.set_zlim3d(5,8)
axagain4.set_ylim3d(1,4)
axagain4.set_xlim3d(0.75,4.5)
axagain4.view_init(elev=35, azim=145)

axagain4.legend()

plt.show()
#********************************************************************************#

aaverage ,baverage,caverage= np.mean(x1),np.mean(y1),np.mean(z1)
daverage ,eaverage,faverage= np.mean(x2),np.mean(y2),np.mean(z2)
gaverage ,haverage,iaverage= np.mean(x3),np.mean(y3),np.mean(z3)



shapesagain7=plt.figure()
axagain7=Axes3D(shapesagain7)
a=axagain7.scatter(aaverage ,baverage,caverage, c='r', label= 'pH 3 Average', marker='D',s=60)
b=axagain7.scatter(daverage ,eaverage,faverage, c='g',label='pH 7 Average',marker='s',s=60)
c=axagain7.scatter(gaverage ,haverage,iaverage, c='purple',label='pH 10 Average',s=60)


axagain7.set_xlabel(deltaaname,fontsize=15) #delatA
axagain7.set_ylabel(normalisedmassname,fontsize=15) #maxmass
axagain7.set_zlabel('$Log_{10}$ Intensity',fontsize=15) #intensity
axagain7.legend()
axagain7.set_zlim3d(5,8)
axagain7.set_ylim3d(1,4)
axagain7.set_xlim3d(0.75,4.5)
axagain7.view_init(elev=35, azim=145)
plt.show()


#********************************************************************************#
# x1 = combodf[ combodf['metal'] == 1 and combodf['time']==3].deltaA.iloc[:]
# y1 = combodf[ combodf['metal'] == 1,combodf['time']==3].meanmass.iloc[:]
# z1 = combodf[ combodf['metal'] == 1,combodf['time']==3].intensity.iloc[:]
# c1 = combodf[ combodf['metal'] == 1,combodf['time']==3].time.iloc[:]

# shapesagain5=plt.figure()
# axagain5=Axes3D(shapesagain5)
# axagain5.scatter(x1,y1,z1,c=c1)
# axagain5.set_xlabel('\u0394A') #delatA
# axagain5.set_ylabel('Mean m/z/ML') #maxmass
# axagain5.set_zlabel('$Log_{10}$ Intensity') #intensity
# # axagain5.set_zlim3d(5,8)
# axagain5.set_ylim3d(1,4)
# axagain5.set_xlim3d(0.75,3.5)


plt.show()
#********************************************************************************#






x1 = combodf[ combodf['cluster']==4].deltaA.iloc[:]
y1 = combodf[ combodf['cluster'] == 4].meanmass.iloc[:]
z1 = combodf[ combodf['cluster'] == 4].intensity.iloc[:]
x2 = combodf[ combodf['cluster']==6].deltaA.iloc[:]
y2 = combodf[ combodf['cluster'] == 6].meanmass.iloc[:]
z2 = combodf[ combodf['cluster'] == 6].intensity.iloc[:]
x3 = combodf[ combodf['cluster'] <= 2].deltaA.iloc[:]
y3 = combodf[ combodf['cluster'] <= 2].meanmass.iloc[:]
z3 = combodf[ combodf['cluster'] <= 2].intensity.iloc[:]

shapesagaing=plt.figure()
axagaing=Axes3D(shapesagaing)
lo=axagaing.scatter(x1,y1,z1, c='orange', label= 'Cluster 1 - Increased Density', marker='D')
mid=axagaing.scatter(x2,y2,z2, c='b',label='Cluster 2 - Increased Mass',marker='s')
hi=axagaing.scatter(x3,y3,z3, c='purple',label='Cluster 3 -  Low Activity')
axagaing.set_xlabel(deltaaname,fontname='Arial',fontsize=25, labelpad=12) #delatA
axagaing.set_ylabel(normalisedmassname,fontname='Arial',fontsize=25, labelpad=12) #maxmass
axagaing.set_zlabel('$Log_{10}$ Intensity',fontname='Arial',fontsize=25, labelpad=12) #intensity
axagaing.set_zlim3d(5,8)
axagaing.set_ylim3d(1,3.5)
axagaing.set_xlim3d(0.75,4.5)
axagaing.zaxis.set_tick_params(labelsize=15)
axagaing.yaxis.set_tick_params(labelsize=15)
axagaing.xaxis.set_tick_params(labelsize=15)
axagaing.view_init(elev=35, azim=145)

axagaing.legend()

plt.show()

#********************************************************************************#

aaverage ,baverage,caverage= np.mean(x1),np.mean(y1),np.mean(z1)
daverage ,eaverage,faverage= np.mean(x2),np.mean(y2),np.mean(z2)
gaverage ,haverage,iaverage= np.mean(x3),np.mean(y3),np.mean(z3)



shapesagainp=plt.figure()
axagainp=Axes3D(shapesagainp)
a=axagainp.scatter(aaverage ,baverage,caverage, c='orange', label= 'CLuster 1 Average', marker='D', s=60)
b=axagainp.scatter(daverage ,eaverage,faverage, c='b',label='Cluster 2 Average',marker='s', s=60)
c=axagainp.scatter(gaverage ,haverage,iaverage, c='purple',label='Cluster 3 Average',s=60)


axagainp.set_xlabel(deltaaname,fontsize=15) #delatA
axagainp.set_ylabel(normalisedmassname,fontsize=15) #maxmass
axagainp.set_zlabel('$Log_{10}$ Intensity',fontsize=15) #intensity
axagainp.legend()
axagainp.set_zlim3d(5,8)
axagainp.set_ylim3d(1,4)
axagainp.set_xlim3d(0.75,4.5)


axagainp.legend()
axagainp.view_init(elev=35, azim=145)
plt.show()

#********************************************************************************#

x1 = combodf[ combodf['Ratio']==3].deltaA.iloc[:]
y1 = combodf[ combodf['Ratio'] == 3].meanmass.iloc[:]
z1 = combodf[ combodf['Ratio'] == 3].intensity.iloc[:]
x2 = combodf[ combodf['Ratio']==2].deltaA.iloc[:]
y2 = combodf[ combodf['Ratio'] == 2].meanmass.iloc[:]
z2 = combodf[ combodf['Ratio'] == 2].intensity.iloc[:]
x3 = combodf[ combodf['Ratio'] <= 1].deltaA.iloc[:]
y3 = combodf[ combodf['Ratio'] <= 1].meanmass.iloc[:]
z3 = combodf[ combodf['Ratio'] <= 1].intensity.iloc[:]
x4 = combodf[ combodf['Ratio'] <= 0.5].deltaA.iloc[:]
y4 = combodf[ combodf['Ratio'] <= 0.5].meanmass.iloc[:]
z4 = combodf[ combodf['Ratio'] <= 0.5].intensity.iloc[:]


shapesagainz=plt.figure()
axagainz=Axes3D(shapesagainz)
lo=axagainz.scatter(x1,y1,z1, c='b', label= 'M:L 3:1', marker='D')
mid=axagainz.scatter(x2,y2,z2, c='yellowgreen',label='M:L 2:1',marker='s')
hi=axagainz.scatter(x3,y3,z3, c='lightcoral',label='M:L 1:1', marker='P')
loow=axagainz.scatter(x4,y4,z4, c='k',label='M"L 0.5:1', marker='h')

axagainz.set_xlabel(deltaaname) #delatA
axagainz.set_ylabel(normalisedmassname) #maxmass
axagainz.set_zlabel('$Log_{10}$ Intensity') #intensity
axagainz.view_init(elev=35, azim=145)

axagaing.set_zlim3d(5,8)
axagaing.set_ylim3d(1,4)
axagaing.set_xlim3d(0.75,4.5)


axagainz.legend()

plt.show()





aaverage ,baverage,caverage= np.mean(x1),np.mean(y1),np.mean(z1)
daverage ,eaverage,faverage= np.mean(x2),np.mean(y2),np.mean(z2)
gaverage ,haverage,iaverage= np.mean(x3),np.mean(y3),np.mean(z3)
javerage ,kaverage,laverage= np.mean(x4),np.mean(y4),np.mean(z4)

shapesagain56=plt.figure(figsize=(10,10))
axagain56=Axes3D(shapesagain56)
a=axagain56.scatter(aaverage ,baverage,caverage, c='b', label= '3:1 Average', marker='D',s=60)
b=axagain56.scatter(daverage ,eaverage,faverage, c='yellowgreen', label= '2:1 Average', marker='s',s=60)
c=axagain56.scatter(gaverage ,haverage,iaverage, c='lightcoral', label= '1:1 Average', marker='P',s=60)
d=axagain56.scatter(javerage ,kaverage,laverage, c='k', label= '0.5:1 Average', marker='h',s=60)


axagain56.set_xlabel(deltaaname,fontsize=15) #delatA
axagain56.set_ylabel(normalisedmassname,fontsize=15) #maxmass
axagain56.set_zlabel('$Log_{10}$ Intensity',fontsize=15) #intensity
axagain56.legend()
axagain56.set_zlim3d(5,8)
axagain56.set_ylim3d(1,4)
axagain56.set_xlim3d(0.75,4.5)

axagain56.legend()
axagain56.view_init(elev=35, azim=145)

plt.show()

#********************************************************************************#

x1 = combodf[ combodf['newclus']==2].deltaA.iloc[:]
y1 = combodf[ combodf['newclus'] == 2].meanmass.iloc[:]
z1 = combodf[ combodf['newclus'] == 2].intensity.iloc[:]
x2 = combodf[ combodf['newclus']==0].deltaA.iloc[:]
y2 = combodf[ combodf['newclus'] == 0].meanmass.iloc[:]
z2 = combodf[ combodf['newclus'] == 0].intensity.iloc[:]
x3 = combodf[ combodf['newclus'] == 1].deltaA.iloc[:]
y3 = combodf[ combodf['newclus'] == 1].meanmass.iloc[:]
z3 = combodf[ combodf['newclus'] == 1].intensity.iloc[:]
# x4 = combodf[ combodf['newclus'] == 3].deltaA.iloc[:]
# y4 = combodf[ combodf['newclus'] == 3].meanmass.iloc[:]
# z4 = combodf[ combodf['newclus'] == 3].intensity.iloc[:]

shapesagainh=plt.figure()
axagainh=Axes3D(shapesagainh)
lo=axagainh.scatter(x1,y1,z1, c='orange', label= 'Cluster 1 - Increased Density', marker='D')
mid=axagainh.scatter(x3,y3,z3, c='b',label='Cluster 2 - Increased Mass',marker='s')
hi=axagainh.scatter(x2,y2,z2, c='purple',label='Cluster 3 - Low Activity')
# extra=axagainh.scatter(x4,y4,z4, c='g',label='Cluster 3 - Low Activity')

axagainh.set_xlabel(deltaaname) #delatA
axagainh.set_ylabel(normalisedmassname) #maxmass
axagainh.set_zlabel('$Log_{10}$ Intensity') #intensity
axagainh.view_init(elev=35, azim=145)

axagainh.set_xlabel(deltaaname,fontname='Arial',fontsize=25, labelpad=12) #delatA
axagainh.set_ylabel(normalisedmassname,fontname='Arial',fontsize=25, labelpad=12) #maxmass
axagainh.set_zlabel('$Log_{10}$ Intensity',fontname='Arial',fontsize=25, labelpad=12) #intensity
axagainh.set_zlim3d(5,8)
axagainh.set_ylim3d(1,3.5)
axagainh.set_xlim3d(0.75,4.5)
axagainh.zaxis.set_tick_params(labelsize=15)
axagainh.yaxis.set_tick_params(labelsize=15)
axagainh.xaxis.set_tick_params(labelsize=15)


axagainh.legend()

plt.show()

x1 = combodf[ combodf['newclus']==2].deltaA.iloc[:]
y1 = combodf[ combodf['newclus'] == 2].meanmass.iloc[:]
z1 = combodf[ combodf['newclus'] == 2].intensity.iloc[:]
x2 = combodf[ combodf['newclus']==0].deltaA.iloc[:]
y2 = combodf[ combodf['newclus'] == 0].meanmass.iloc[:]
z2 = combodf[ combodf['newclus'] == 0].intensity.iloc[:]
x3 = combodf[ combodf['newclus'] == 1].deltaA.iloc[:]
y3 = combodf[ combodf['newclus'] == 1].meanmass.iloc[:]
z3 = combodf[ combodf['newclus'] == 1].intensity.iloc[:]

shapesagaino=plt.figure()
axagaino=Axes3D(shapesagaino)
# lo=axagainh.scatter(x1,y1,z1, c='r', label= 'Cluster 1- Increased Density', marker='D')
# mid=axagainh.scatter(x2,y2,z2, c='g',label='Cluster 2 - Lower Activity',marker='s')
# hi=axagainh.scatter(x3,y3,z3, c='purple',label='Cluster 3-  High Linearity')

axagaino.set_xlabel(deltaaname) #delatA
axagaino.set_ylabel(normalisedmassname) #maxmass
axagaino.set_zlabel('$Log_{10}$ Intensity') #intensity
axagaino.view_init(elev=35, azim=145)

axagaing.set_zlim3d(5,8)
axagaing.set_ylim3d(1,4)
axagaing.set_xlim3d(0.75,4.5)


axagainh.legend()

plt.show()


x1 = combodf[ combodf['time']==0].deltaA.iloc[:]
y1 = combodf[ combodf['time'] == 0].meanmass.iloc[:]
z1 = combodf[ combodf['time'] == 0].intensity.iloc[:]
x2 = combodf[ combodf['time']==7].deltaA.iloc[:]
y2 = combodf[ combodf['time'] == 7].meanmass.iloc[:]
z2 = combodf[ combodf['time'] == 7].intensity.iloc[:]
x3 = combodf[ combodf['time'] == 100].deltaA.iloc[:]
y3 = combodf[ combodf['time'] == 100].meanmass.iloc[:]
z3 = combodf[ combodf['time'] == 100].intensity.iloc[:]

shapesagainl=plt.figure()
axagainl=Axes3D(shapesagainl)
lo=axagainl.scatter(x1,y1,z1, c='r', label= 'Start', marker='D')
mid=axagainl.scatter(x2,y2,z2, c='g',label='Mid',marker='s')
hi=axagainl.scatter(x3,y3,z3, c='purple',label='End')





axagainl.set_xlabel(deltaaname) #delatA
axagainl.set_ylabel(normalisedmassname) #maxmass
axagainl.set_zlabel('$Log_{10}$ Intensity') #intensity
axagainl.view_init(elev=35, azim=145)

axagainl.set_zlim3d(5,8)
axagainl.set_ylim3d(1,4)
axagainl.set_xlim3d(0.75,4.5)


axagainl.legend()

plt.show()


x1 = combodf[ combodf['Ligandno'] == 1].deltaA.iloc[:]
y1 = combodf[ combodf['Ligandno'] == 1].meanmass.iloc[:]
z1 = combodf[ combodf['Ligandno'] == 1].intensity.iloc[:]
x2 = combodf[ combodf['Ligandno'] == 2].deltaA.iloc[:]
y2 = combodf[ combodf['Ligandno'] == 2].meanmass.iloc[:]
z2 = combodf[ combodf['Ligandno'] == 2].intensity.iloc[:]
x3 = combodf[ combodf['Ligandno'] == 3].deltaA.iloc[:]
y3 = combodf[ combodf['Ligandno'] == 3].meanmass.iloc[:]
z3 = combodf[ combodf['Ligandno'] == 3].intensity.iloc[:]
x4 = combodf[ combodf['Ligandno'] == 4].deltaA.iloc[:]
y4 = combodf[ combodf['Ligandno'] == 4].meanmass.iloc[:]
z4 = combodf[ combodf['Ligandno'] == 4].intensity.iloc[:]
x5 = combodf[ combodf['Ligandno'] == 5].deltaA.iloc[:]
y5 = combodf[ combodf['Ligandno'] == 5].meanmass.iloc[:]
z5 = combodf[ combodf['Ligandno'] == 5].intensity.iloc[:]
x6 = combodf[ combodf['Ligandno'] == 6].deltaA.iloc[:]
y6 = combodf[ combodf['Ligandno'] == 6].meanmass.iloc[:]
z6 = combodf[ combodf['Ligandno'] == 6].intensity.iloc[:]



shapesagainpp=plt.figure()
axagainpp=Axes3D(shapesagainpp)
c=axagainpp.scatter(x3,y3,z3, c='purple',label='N/A')
a=axagainpp.scatter(x1,y1,z1, c='y', label= '180', marker='D')
d=axagainpp.scatter(x4,y4,z4, c='y',label='180',marker='D')
b=axagainpp.scatter(x2,y2,z2, c='k',label='60',marker='s')
f=axagainpp.scatter(x6,y6,z6, c='b',label='0', marker='P')
e=axagainpp.scatter(x5,y5,z5, c='b',label='0',marker="P")

axagainpp.set_xlabel(deltaaname) #delatA
axagainpp.set_ylabel(normalisedmassname) #maxmass
axagainpp.set_zlabel('$Log_{10}$ Intensity') #intensity
axagainpp.set_zlim3d(5,8)
axagainpp.set_ylim3d(1,4)
axagainpp.set_xlim3d(0.75,3.5)

handles, labels = axagain2.get_legend_handles_labels()
totallen=len(x1)+len(x2)+len(x3)+len(x4)+len(x5)+len(x6)
patch = mpatches.Patch(color='k', label=(str(totallen)+' Spectra'))
handles.append(patch)
axagainpp.legend(handles=handles)

axagainpp.view_init(elev=35, azim=145)
plt.show()



# x1 = combodf[ combodf['metal'] == 2].deltaA.iloc[:]
# y1 = combodf[ combodf['metal'] == 2].meanmass.iloc[:]
# z1 = combodf[ combodf['metal'] == 2].intensity.iloc[:]




# shapesagainvvv=plt.figure()
# axagainvvv=Axes3D(shapesagainvvv)

# a=axagainvvv.scatter(x1,y1,z1, c=combodf[ combodf['metal'] == 2].expno.iloc[:], marker='D')

# axagainvvv.set_xlabel(deltaaname) #delatA
# axagainvvv.set_ylabel(normalisedmassname) #maxmass
# axagainvvv.set_zlabel('$Log_{10}$ Intensity') #intensity
# axagainvvv.set_zlim3d(5,8)
# axagainvvv.set_ylim3d(1,4)
# axagainvvv.set_xlim3d(0.75,3.5)

# handles, labels = axagainvvv.get_legend_handles_labels()
# totallen=len(x1)
# patch = mpatches.Patch(color='k', label=(str(totallen)+' Spectra'))
# handles.append(patch)
# axagainvvv.legend(handles=handles)

# axagainvvv.view_init(elev=35, azim=145)
# plt.show()


# for value in range(72):
#     REP=plt.figure()
#     rep=Axes3D(REP)
#     x1 = combodf[ (combodf['metal'] == 0) &( combodf['expno']==value) ].deltaA.iloc[:]
#     y1 = combodf[ (combodf['metal'] == 0) &( combodf['expno']==value) ].meanmass.iloc[:]
#     z1 = combodf[ (combodf['metal'] == 0) &( combodf['expno']==value) ].intensity.iloc[:]
#     x2 = combodf[ (combodf['metal'] == 1) &( combodf['expno']==value) ].deltaA.iloc[:]
#     y2 = combodf[ (combodf['metal'] == 1) &( combodf['expno']==value) ].meanmass.iloc[:]
#     z2 = combodf[ (combodf['metal'] == 1) &( combodf['expno']==value) ].intensity.iloc[:]
#     x3 = combodf[ (combodf['metal'] == 2) &( combodf['expno']==value) ].deltaA.iloc[:]
#     y3 = combodf[ (combodf['metal'] == 2) &( combodf['expno']==value) ].meanmass.iloc[:]
#     z3 = combodf[ (combodf['metal'] == 2) &( combodf['expno']==value) ].intensity.iloc[:]
    
    
#     c=rep.scatter(x3,y3,z3, c=combodf[ (combodf['metal'] == 0) &( combodf['expno']==value) ].time.iloc[:],label='Iron')
#     a=rep.scatter(x1,y1,z1, c=combodf[ (combodf['metal'] == 0) &( combodf['expno']==value) ].time.iloc[:], label= 'Copper', marker='D')
#     b=rep.scatter(x2,y2,z2, c=combodf[ (combodf['metal'] == 0) &( combodf['expno']==value) ].time.iloc[:],label='Nickel',marker='s')
    
    
#     rep.set_xlabel(deltaaname) #delatA
#     rep.set_ylabel(normalisedmassname) #maxmass
#     rep.set_zlabel('$Log_{10}$ Intensity') #intensity
#     rep.set_zlim3d(5,8)
#     rep.set_ylim3d(1,4)
#     rep.set_xlim3d(0.75,3.5)
    
#     handles, labels = rep.get_legend_handles_labels()
#     totallen=len(x1)+len(x2)+len(x3)
#     patch = mpatches.Patch(color='k', label=(str(totallen)+' Spectra'))
#     handles.append(patch)
#     rep.legend(handles=handles)
    
#     rep.view_init(elev=35, azim=145)
#     plt.show()


REP=plt.figure()
rep=Axes3D(REP)
x1 = combodf[ ( combodf['Ligandno']==1) |( combodf['Ligandno']==4) ].deltaA.iloc[:]
y1 = combodf[ ( combodf['Ligandno']==1) |( combodf['Ligandno']==4) ].meanmass.iloc[:]
z1 = combodf[ ( combodf['Ligandno']==1) |( combodf['Ligandno']==4) ].intensity.iloc[:]
x2 = combodf[ ( combodf['Ligandno']==2) |( combodf['Ligandno']==6) ].deltaA.iloc[:]
y2 = combodf[ ( combodf['Ligandno']==2) |( combodf['Ligandno']==6) ].meanmass.iloc[:]
z2 = combodf[ ( combodf['Ligandno']==2) |( combodf['Ligandno']==6) ].intensity.iloc[:]
x3 = combodf[ ( combodf['Ligandno']==5)  ].deltaA.iloc[:]
y3 = combodf[ ( combodf['Ligandno']==5)  ].meanmass.iloc[:]
z3 = combodf[ ( combodf['Ligandno']==5)  ].intensity.iloc[:]
x4 = combodf[ ( combodf['Ligandno']==3)  ].deltaA.iloc[:]
y4 = combodf[ ( combodf['Ligandno']==3)  ].meanmass.iloc[:]
z4 = combodf[ ( combodf['Ligandno']==3)  ].intensity.iloc[:]
    
    
c=rep.scatter(x3,y3,z3, c='y',label='0')
a=rep.scatter(x1,y1,z1, c='b', label= '180', marker='D')
b=rep.scatter(x2,y2,z2, c='g',label='60',marker='s')
d=rep.scatter(x4,y4,z4, c='m',label='None',marker='s')   
    
rep.set_xlabel(deltaaname) #delatA
rep.set_ylabel(normalisedmassname) #maxmass
rep.set_zlabel('$Log_{10}$ Intensity') #intensity
rep.set_zlim3d(5,8)
rep.set_ylim3d(1,4)
rep.set_xlim3d(0.75,3.5)
    
handles, labels = rep.get_legend_handles_labels()
totallen=len(x1)+len(x2)+len(x3)
patch = mpatches.Patch(color='k', label=(str(totallen)+' Spectra'))
handles.append(patch)
rep.legend(handles=handles)
    
rep.view_init(elev=35, azim=145)
plt.show()  
