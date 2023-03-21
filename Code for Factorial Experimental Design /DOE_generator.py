# -*- coding: utf-8 -*-
"""
Spyder Editor
"""



##This file was written by OLW to generate experimental design for "Automated Structural Activity Screening of β-Diketonate Assemblies with High-Throughput Ion Mobility - Mass Spectrometry"


import numpy as np
import pyDOE2
import pandas as pd

##Factors

#Metal Choice here
#in this case Copper acac 
metals=['Cuacac']

#ligands

ligands=['dahd','dhba','f6acac',  'DHNQ','ddpdp','ppd']
metal_ratio=[3,2,1,0.5]

#pH
#3 pH points chosen
pH=['Low', 'High','Unaltered']

#volume. End experiment volume will be 1 mL but some pH modification will need volume.
volume=[950]

#how many experiments are there
levels=[len(metals),len(ligands),len(metal_ratio),len(pH),len(volume)]

#generate experiment using number of variables contained in levels
DOE=pyDOE2.fullfact(levels)
DOEdf=pd.DataFrame(DOE,columns=['Metal','Ligand','Metal Ratio','pH','Volume']) #heating removed

#renaming generated DOE datafram values with text names for readability
DOEdf["Metal"].replace({0: metals[0]}, inplace=True)
DOEdf["Ligand"].replace({0: ligands[0], 1:ligands[1],2:ligands[2],3:ligands[3],4:ligands[4],5:ligands[5]}, inplace=True)
DOEdf["Metal Ratio"].replace({0: metal_ratio[0], 1:metal_ratio[1],2:metal_ratio[2],3:metal_ratio[3]}, inplace=True)
DOEdf["pH"].replace({0: pH[0],1: pH[1],2: pH[2]}, inplace=True)
DOEdf["Volume"].replace({0: volume[0]}, inplace=True)

#printing and saving      
print(DOEdf)
DOEdf.to_csv('CuDOE.csv')