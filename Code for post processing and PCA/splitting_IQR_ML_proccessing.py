# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 12:03:27 2021

@author: oscar
"""

"""
Created on Tue May 18 14:26:p29 2021

@author: oscar
"""

import numpy as np

import pandas as pd
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D

import scipy
#from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.cluster import AgglomerativeClustering
# from sklearn.linear_model import LinearRegression
# from sklearn.cluster import DBSCAN
# from sklearn.cluster import OPTICS, cluster_optics_dbscan
# import statsmodels.api as sm

# In[54]:

    
plt.rc('legend',fontsize=15)
font={'size':15,}
plt.rc('font',family='Arial')
plt.rcParams['axes.linewidth'] = 1

##*********** Read in the DOE file
csv_ids = pd.read_csv('C:/Users/oscar/OneDrive - UNSW/Paperdraft/DOE.csv')


# In[55]:


##***********Read in the xls
xl = pd.ExcelFile('C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/Cu_D10_trunc.xls')
outname="C:/Users/oscar/OneDrive - UNSW/Paperdraft/PaperdraftMS/IQRd/outtest.csv"
time=10
metalinp="Ni"

# print(xl.sheet_names)
# then you can assign them to a dict



# d = {} # your dict.
all_data = {}
for sheet in xl.sheet_names:
    #print(sheet)
    all_data[f'{sheet}']= pd.read_excel(xl,sheet_name=sheet)
#     d[f'{sheet}']= pd.read_excel(xl,sheet_name=sheet)
#print(all_data)
# print(d.values())
###**********This is used to clean-up the data
samples_ids = list(all_data.keys())

metals={'Cu':62.93,'Ni':57.935347,'Fe':55.845, 'Zn':65.38}
ligands={"ddpdp":328.17,"dahd":196.07,"dhba":152.01, "DHNQ":188.01, "f6acac":206.99, "ppd":244.07}

frog=[]
names=[]

for value in all_data.values():
    df=value
    del df['Observed mass (Da)']
    del df['Expected mass (Da)']
    del df['Mass error (mDa)']
    del df['Mass error (ppm)']
    del df['Expected RT (min)']
    del df['Observed RT (min)']
   
    del df['Expected CCS (Ǻ2)']
    del df['CCS delta (%)']
    
    del df['Adducts']
    
    

    
    df['dtwieghted']=df['Observed drift (ms)']**0.894
    
ymean_all=[]
rmse_all=[]

output=[]

for ids in samples_ids:
        
        extract = all_data[ids]
        #extract['Multiply']=extract['Observed m/z']*extract['Response']
        sample_number = ids.split('_')
        if sample_number[0] !="Blank":
            if sample_number[0] !="X":
            
            
                ssss = sample_number[1]
                lookligand = np.asarray(csv_ids[csv_ids['sampleno'] == int(ssss)]['Ligand'])
                lookpH = np.asarray(csv_ids[csv_ids['sampleno'] == int(ssss)]['pH'])
                pH=lookpH[0]
                ratiolookup=np.asarray(csv_ids[csv_ids['sampleno'] == int(ssss)]['Metal Ratio'])
                ratio=ratiolookup[0]
                
                ligandm=ligands.get(lookligand[0])
                metalm=float(metals.get(metalinp))
                ml=float(metalm)+float(ligandm)
                
                extract['m/z/ml']=extract['Observed m/z']/ml
                #extract['caledCCS']=(extract['Observed CCS (Ǻ2)'])*1.0301+1.2335
                intens=extract['Detector counts'].sum(skipna=True)
                # extract = extract.drop(extract[extract['Detector counts'] < (intens/1000)].index, inplace = True)
                dfz = extract[extract['Detector counts'] > (intens/10000)]  
                
                meanmass=extract['m/z/ml'].mean()
                medianmass=extract['m/z/ml'].median()
                maxmass=extract['m/z/ml'].max()
                
                # plt.figure()
                plt.subplots(1,2,figsize=(10,4))
                
                
                
                
                data=[extract['m/z/ml'],extract['dtwieghted']]
                headers=['m/z/ml','dtwieghted']
                   
                
                
                df3 = pd.concat(data, axis=1, keys=headers)
                df4= df3.dropna(how='any')
                 
                
                X=df4['m/z/ml']
                
                Y=df4['dtwieghted']
                
               
    
                
                #plt.suptitle(str(ids)+" "+str(lookligand))
                
                # clusno=len(df4["Cluster"].value_counts())
                # l = [i for i in range(clusno)]
                
                
                
                colorkey={0:'orange', 1:'green',2:'black', 3:'blue',4:"yellow"}
                
                plt.subplot(121)
                plt.ylabel('Drift Time$^{0.894}$')
                plt.xlabel('$\it{m/z}$/ml')
                
                            #df4['Cluster'])
                plt.ylim(0,15)
                plt.xlim(0,8)
                
                
                
               
                
                X=X.to_numpy()
                Y=Y.to_numpy()
                
                dcX=X
                dcY=Y
                
                first_quartile = np.quantile(dcY, 0.25)
                third_quartile = np.quantile(dcY, 0.75)
                
                ind_dcYlows=[]
                ind_dcYhighs=[]
                
                lowouts=first_quartile-1.5*(third_quartile-first_quartile)
                highouts=third_quartile+1.5*(third_quartile-first_quartile)
                
                ind_dcYlows=np.where(dcY<lowouts)
                ind_dcYhighs=np.where(dcY>highouts)
                ind_dcYrem=np.concatenate((ind_dcYlows,ind_dcYhighs),axis=1)
                
                removedX=dcX[ind_dcYrem]
                removedY=dcY[ind_dcYrem]
                
                
                
                dcX=np.delete(dcX,ind_dcYrem)
                dcY=np.delete(dcY,ind_dcYrem)

                
                a1,b1=np.polyfit(dcX,dcY,1)
                
                dcypred=a1*dcX+b1
                
                
                plt.subplot(121)
                maxval=math.ceil(np.max(dcX))
                z = np.linspace(0, maxval)
                #plt.plot(z, b1+(a1*z),'k')
                
                a1round=round(a1,3)
                
                
                dcerror=dcY-dcypred
                rmse_tmp=np.sqrt(np.mean(np.power(dcerror,2)))
                rmse_all.append(rmse_tmp)
                
                ymean_all.append(np.mean(dcY))
                
                
                
                plt.subplot(122)
                
                a=plt.hist(dcerror,100, density=1, color='grey')
                plt.ylabel('Error Density')
                plt.xlabel('Root Mean Squared Error')
                
                mu, sigma = scipy.stats.norm.fit(dcerror)
                best_fit_line = scipy.stats.norm.pdf(a[1], mu, sigma)
                bin_width = a[1][1] - a[1][0]
                area = np.sum(a[0]) * bin_width
                
                bin_mids = (a[1][:-1] + a[1][1:]) / 2
                squared_error = ((scipy.stats.norm.pdf(bin_mids, mu, sigma) * area - a[0]) ** 2).sum()
                #print("squared_error div degrees_of_freedom: ", squared_error / (a[0].size - 3))
                
                
                #plt.plot(a[1], best_fit_line,"-.k")
                #plt.title(squared_error)
                
                plt.xlim(-4,4)
                
                
                
                
                
                
                if rmse_tmp/np.mean(dcY) > 0.100:  #0.125
                
                    
                    
                    
                    #skimming opp pecentiles from histogram to get narrow version for local minimum
                    hist_dcerror=np.concatenate((a[0][:,None],a[1][:-1,None]),axis=1)
                    hist_dcerror_10_90=hist_dcerror[(hist_dcerror[:,1]>np.percentile(dcerror,10)) & (hist_dcerror[:,1]<np.percentile(dcerror,75)), :]
                    ind_min_hist_10_90=np.nonzero(hist_dcerror_10_90[:,0] == np.min(hist_dcerror_10_90[:,0]))
                    
                    
                    
                    #splitting at minima twin peaks
                    
                    ind_min_peak= np.nonzero(dcerror<(hist_dcerror_10_90[ind_min_hist_10_90[0],1][0]))
                    ind_max_peak= np.nonzero(dcerror>=(hist_dcerror_10_90[ind_min_hist_10_90[0],1][0]))
                    
                    a_min_peak,b_min_peak=np.polyfit(dcX[ind_min_peak],dcY[ind_min_peak],1)
                    pred_min_peak=a_min_peak*(dcX[ind_min_peak])+b_min_peak
                    
                    a_max_peak,b_max_peak=np.polyfit(dcX[ind_max_peak],dcY[ind_max_peak],1)
                    pred_max_peak=a_max_peak*(dcX[ind_max_peak])+b_max_peak
                    
                    
                    
                    
                    
                    
                    
                    plt.xlim(-4,4)
                    
                    
                    
                    #print("ID:"+str(ids)+" A "+str(a1round))
                    
                    
                    plt.subplot(121)
                    z = np.linspace(0, maxval)
                    plt.plot(z, b1+(a1*z),'--k')
                    
                    
                    plt.plot(z, b_max_peak+a_max_peak*(z),c='k')
                    plt.scatter(dcX[ind_max_peak],dcY[ind_max_peak],c='orange',marker="^")
                    
                    plt.plot(z, b_min_peak+a_min_peak*(z),c='k')
                    plt.scatter(dcX[ind_min_peak],dcY[ind_min_peak],c='blue',marker="v")
                    
                    plt.scatter(removedX,removedY,c='grey',alpha=0.7)
                    
                    dcypred_dual=np.empty([len(dcY),])
                    
                    dcypred_dual[:,]=np.nan
                    dcypred_dual[ind_max_peak]=pred_max_peak
                    dcypred_dual[ind_min_peak]=pred_min_peak                    
                    
                    
                        
                    newerror=dcY-dcypred_dual
                    plt.subplot(122)
                    
                    
                    b=plt.hist(newerror,100,density=1,alpha=0.7,color='green')
                    # ce=plt.hist(dcerror[ind_max_peak],100,density=1,alpha=1,color='orange')
                    # de=plt.hist(dcerror[ind_min_peak],100,density=1,alpha=1,color='b')
                    
                    newmu, newsigma = scipy.stats.norm.fit(newerror)
                    best_fit_line_new = scipy.stats.norm.pdf(a[1], newmu, newsigma,)
                    #plt.plot(a[1], best_fit_line_new,'k')
                    
                    #plt.plot([np.percentile(dcerror,10),np.percentile(dcerror,10)],[0,np.max(b[0])],'--k')
                    
                    #plt.plot([np.percentile(dcerror,75),np.percentile(dcerror,75)],[0,np.max(b[0])],'-.k')
                    
                    
                    
                    a_max_peak_round=round(a_max_peak,3)
                    a_min_peak_round=round(a_min_peak,3)
                    
                    outputlist=(ids, "Split", a1round,a_max_peak_round, a_min_peak,lookligand[0],metalinp,12,np.log10(intens),meanmass,medianmass,maxmass,pH,lookligand[0],ratio)
                
                                    
                
                
                
                
                else:
                    
                    #plt.plot([np.percentile(dcerror,10),np.percentile(dcerror,10)],[0,np.max(a[0])],'--k')
                #plt.plot([np.percentile(dcerror,90),np.percentile(dcerror,90)],[0,np.max(a[0])],'k')
                    
                #plt.plot([np.percentile(dcerror,25),np.percentile(dcerror,25)],[0,np.max(a[0])],'--k')
                    #plt.plot([np.percentile(dcerror,75),np.percentile(dcerror,75)],[0,np.max(a[0])],'-.k')
                    
                    plt.subplot(121)
                    #z = np.linspace(0, 7)
                    
                    plt.scatter(x=X,y=(Y),c='grey')
                    plt.plot(z, b1+(a1*z),'--r')
                    
                    outputlist=(ids, "NoSplit", a1round,a1round,a1round,lookligand[0],metalinp,12,np.log10(intens),meanmass,medianmass,maxmass,pH,lookligand[0],ratio)
                
                output.append(outputlist)
                
                plt.show()
                
outdf=pd.DataFrame(output)
outdf=outdf.rename({0:"sample",1:"IsSplit?",2:"Aorg",3:"Aabove",4:"Abelow",5:"Ligand",6:"metal",8:'intensity',9:'meanmass',10:'medianmass',11:"maxmass",12:"pH",13:"Ligandno",14:"Ratio"},axis="columns")
                
outdf['deltaA']=outdf['Aabove'].dropna()/outdf['Abelow'].dropna()
outdf['Ligandno'] = outdf['Ligandno'].replace(['dahd', 'dhba' ,'f6acac' ,'DHNQ' ,'ddpdp' ,'ppd'],[1,2,3,4,5,6])                 
outdfsave=outdf.copy()


shapes2=plt.figure()

ax2=Axes3D(shapes2)
ax2.scatter(outdf['maxmass'],outdf['meanmass'],outdf['medianmass'])
ax2.set_xlabel('maxmass')
ax2.set_ylabel('meanmass')
ax2.set_zlabel('medianmass')


plt.show()

newdf=outdf.copy()


del newdf['sample']
del newdf['metal']
del newdf['medianmass']
del newdf['maxmass']
del newdf['Aabove']
del newdf['Abelow']
del newdf['Aorg']
del newdf['pH']
del newdf['Ligand']
del newdf['Ligandno']
del newdf['Ratio']
          
newdf['IsSplit?'] = newdf['IsSplit?'].replace(['NoSplit','Split'],[0,1])            
outdf['pH'] = outdf['pH'].replace(['High','Unaltered',"Low"],[10,7,3]) 

              
cluster = AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
cluster.fit_predict(newdf)


                
outdf['cluster']=cluster.labels_
outdf['time']=time

shapes=plt.figure()

ax=Axes3D(shapes)

xs=newdf['deltaA']
ys=newdf['meanmass']
zs=newdf['intensity']


ax.scatter(xs,ys,zs,c=cluster.labels_)
ax.set_xlabel('deltaA') #delatA
ax.set_ylabel('meanmass') #maxmass
ax.set_zlabel('intensity') #intensity

# for i in range(outdf.shape[0]):
#       x = xs[i]
#       y = ys[i]
#       z = zs[i]
#       label = i
#       ax.scatter(x, y, z, color='b')
#       ax.text(x, y, z, '%s' % (label), size=10, zorder=1, color='k')

plt.show()

shapesagain=plt.figure()

xs=newdf['deltaA']
ys=newdf['meanmass']
zs=newdf['intensity']
cs=outdf['pH']

axagain=Axes3D(shapesagain)
axagain.scatter(xs,ys,zs,c=cs)
axagain.set_xlabel('deltaA') #delatA
axagain.set_ylabel('meanmass') #maxmass
axagain.set_zlabel('intensity') #intensity

# for i in range(outdf.shape[0]):
#       x = xs[i]
#       y = ys[i]
#       z = zs[i]
#       label = i
#       ax.scatter(x, y, z, color='b')
#       ax.text(x, y, z, '%s' % (label), size=10, zorder=1, color='k')

plt.show()

xs=newdf['deltaA']
ys=newdf['meanmass']
zs=newdf['intensity']
cs=outdf['Ligandno']
shapesagain2=plt.figure()
axagain2=Axes3D(shapesagain2)
axagain2.scatter(xs,ys,zs,c=cs)
axagain2.set_xlabel('deltaA') #delatA
axagain2.set_ylabel('meanmass') #maxmass
axagain2.set_zlabel('intensity') #intensity

# for i in range(outdf.shape[0]):
#       x = xs[i]
#       y = ys[i]
#       z = zs[i]
#       label = i
#       ax.scatter(x, y, z, color='b')
#       ax.text(x, y, z, '%s' % (label), size=10, zorder=1, color='k')

plt.show()

#outdf['Ligand'] = outdf['Ligand'].replace([1,2,3,4,5,6],['dahd', 'dhba' ,'f6acac' ,'DHNQ' ,'ddpdp' ,'ppd'])   
outdf.to_csv(outname)



                    
                