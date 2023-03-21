# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 15:18:06 2020

@author: oscar
"""

## This file was written by OLW tocreate the experiments designed in DOE_generator.py
## This is for the manuscript "Automated Structural Activity Screening of β-Diketonate Assemblies with High-Throughput Ion Mobility - Mass Spectrometry"
## The OT-2 liquid handling robot was used to prepare these samples.
## Defintions for custom labware are provided.
## Pippettes used were 1000 uL and 300 uL
##metal in this case is Nickel

from opentrons import protocol_api
metadata = {'apiLevel':'2.7',
           'author':'Oscar Lloyd Williams'}

import csv
import subprocess


## text copied across, rather than uploading CSV to OT-2.
csv_copy= '''
0,Niacac,dahd,3.0,Low,950.0
1,Niacac,dhba,3.0,Low,950.0
2,Niacac,f6acac,3.0,Low,950.0
3,Niacac,DHNQ,3.0,Low,950.0
4,Niacac,ddpdp,3.0,Low,950.0
5,Niacac,ppd,3.0,Low,950.0
6,Niacac,dahd,2.0,Low,950.0
7,Niacac,dhba,2.0,Low,950.0
8,Niacac,f6acac,2.0,Low,950.0
9,Niacac,DHNQ,2.0,Low,950.0
10,Niacac,ddpdp,2.0,Low,950.0
11,Niacac,ppd,2.0,Low,950.0
12,Niacac,dahd,1.0,Low,950.0
13,Niacac,dhba,1.0,Low,950.0
14,Niacac,f6acac,1.0,Low,950.0
15,Niacac,DHNQ,1.0,Low,950.0
16,Niacac,ddpdp,1.0,Low,950.0
17,Niacac,ppd,1.0,Low,950.0
18,Niacac,dahd,0.5,Low,950.0
19,Niacac,dhba,0.5,Low,950.0
20,Niacac,f6acac,0.5,Low,950.0
21,Niacac,DHNQ,0.5,Low,950.0
22,Niacac,ddpdp,0.5,Low,950.0
23,Niacac,ppd,0.5,Low,950.0
24,Niacac,dahd,3.0,High,950.0
25,Niacac,dhba,3.0,High,950.0
26,Niacac,f6acac,3.0,High,950.0
27,Niacac,DHNQ,3.0,High,950.0
28,Niacac,ddpdp,3.0,High,950.0
29,Niacac,ppd,3.0,High,950.0
30,Niacac,dahd,2.0,High,950.0
31,Niacac,dhba,2.0,High,950.0
32,Niacac,f6acac,2.0,High,950.0
33,Niacac,DHNQ,2.0,High,950.0
34,Niacac,ddpdp,2.0,High,950.0
35,Niacac,ppd,2.0,High,950.0
36,Niacac,dahd,1.0,High,950.0
37,Niacac,dhba,1.0,High,950.0
38,Niacac,f6acac,1.0,High,950.0
39,Niacac,DHNQ,1.0,High,950.0
40,Niacac,ddpdp,1.0,High,950.0
41,Niacac,ppd,1.0,High,950.0
42,Niacac,dahd,0.5,High,950.0
43,Niacac,dhba,0.5,High,950.0
44,Niacac,f6acac,0.5,High,950.0
45,Niacac,DHNQ,0.5,High,950.0
46,Niacac,ddpdp,0.5,High,950.0
47,Niacac,ppd,0.5,High,950.0
48,Niacac,dahd,3.0,Unaltered,1000.0
49,Niacac,dhba,3.0,Unaltered,1000.0
50,Niacac,f6acac,3.0,Unaltered,1000.0
51,Niacac,DHNQ,3.0,Unaltered,1000.0
52,Niacac,ddpdp,3.0,Unaltered,1000.0
53,Niacac,ppd,3.0,Unaltered,1000.0
54,Niacac,dahd,2.0,Unaltered,1000.0
55,Niacac,dhba,2.0,Unaltered,1000.0
56,Niacac,f6acac,2.0,Unaltered,1000.0
57,Niacac,DHNQ,2.0,Unaltered,1000.0
58,Niacac,ddpdp,2.0,Unaltered,1000.0
59,Niacac,ppd,2.0,Unaltered,1000.0
60,Niacac,dahd,1.0,Unaltered,1000.0
61,Niacac,dhba,1.0,Unaltered,1000.0
62,Niacac,f6acac,1.0,Unaltered,1000.0
63,Niacac,DHNQ,1.0,Unaltered,1000.0
64,Niacac,ddpdp,1.0,Unaltered,1000.0
65,Niacac,ppd,1.0,Unaltered,1000.0
66,Niacac,dahd,0.5,Unaltered,1000.0
67,Niacac,dhba,0.5,Unaltered,1000.0
68,Niacac,f6acac,0.5,Unaltered,1000.0
69,Niacac,DHNQ,0.5,Unaltered,1000.0
70,Niacac,ddpdp,0.5,Unaltered,1000.0
71,Niacac,ppd,0.5,Unaltered,1000.0
'''


#reading csv data into a list

csv_data = csv_copy.splitlines()[1:] # Discard the blank first line.
csv_reader = csv.reader(csv_data)
csv_list=list(csv_reader)


## all OT-2 runs run as functions
def run(protocol: protocol_api.ProtocolContext):
        
    #Defining the location of the methanol bottles
    solventsource=protocol.load_labware('olw_2schott_100000ul',2)
    solvents={'MeOH':'A1', 'MeOH2':'A2'}
    
    #Defining the HPLC vial containing plates for final solutions
    finplate1 = protocol.load_labware('olw_40hplc_2000ul', 4)
    finplate2 = protocol.load_labware('olw_40hplc_2000ul', 5)
    
    #Defining the scintillation vial containing plates for stock solutions
    stockplate=protocol.load_labware('olw_8_wellplate_20000ul', 1)
    stocks={'Niacac':'A1','Niacac2':'A2', 'dahd':'A3','dhba':'A4','f6acac':'B1','DHNQ':'B2','ddpdp':'B3','ppd':'B4'}
    
    
    #Defining locations for tipracks for the 2 pippettes
    tiprackbig = protocol.load_labware('opentrons_96_tiprack_1000ul', 3)
    tipracklil= protocol.load_labware('opentrons_96_tiprack_300ul', 6)
    
    #Defining pippettes and modifying clearence from bottom to avoid errors with inconsistently printed labware.
    pipettebig = protocol.load_instrument('p1000_single_gen2', 'left',tip_racks=[tiprackbig])
    pipettesmol = protocol.load_instrument('p300_single_gen2', 'right',tip_racks=[tipracklil])
    pipettebig.well_bottom_clearance.dispense = 3
    
    #Performing a test transfer to verify protocol setup is correct.
    pipettesmol.pick_up_tip()
    pipettesmol.aspirate(100,stockplate['A1'])
    pipettesmol.dispense(100,finplate2.wells()[39])
    pipettesmol.drop_tip()
    pipettebig.pick_up_tip()
    pipettebig.aspirate(900,solventsource['A1'])
    pipettebig.dispense(900,finplate2.wells()[39])
    pipettebig.drop_tip()
    protocol.pause('Check for success')
    
    #all rows get a constant amount of metal added first
    
    pipettesmol.pick_up_tip()
    
    for row in csv_list:
            finplate=finplate1
            stockposition=str(stocks['Niacac'])
            position=float(row[0])
            if (position  % 2) ==0:
                stockposition=str(stocks['Niacac2'])
            if position >79:
                position=position-80
               
            if position >39:
                position=position-40
                finplate=finplate2
            fevolume=int(100)
            position2=int(position)
            pipettesmol.aspirate(fevolume,stockplate[stockposition])
            pipettesmol.dispense(fevolume,finplate.wells()[position2])
            row[5]=float(row[5])-fevolume
                
    #drop iron tip       
    pipettesmol.drop_tip()
    protocol.comment('Iron Done')
    
    #for loop to deposit ligand stock. Volume is calculated off 100 uL x ratio. 
    for row in csv_list:
        finplate=finplate1
        ligand=row[2]
        pipettesmol.pick_up_tip()
        source=stockplate
        stockposition=str(stocks[ligand])
        position=float(row[0])
        
                
        if position >39:
            position=position-40
            finplate=finplate2
        row3=float(row[3])
        position2=int(position)
        ligvolume=(100*row3)
        pipettesmol.aspirate(ligvolume,source[stockposition])
        pipettesmol.dispense(ligvolume,finplate.wells()[position2])
        pipettesmol.drop_tip()
        row[5]=float(row[5])-ligvolume
        
    protocol.comment('Ligands Done')      
     
    
    # for loop to deposit methanol. Volume calculated by subtracting total amount already pippetted from final amount as given in DOE. 
    # some final volume values are 950 and some are 1000. This is because some were manually modified for pH.
    for row in csv_list:
        finplate=finplate1
        pipettebig.pick_up_tip()
        position= int(row[0])
       
        if position >39:
            position=position-40
            finplate=finplate2
        
        mvolume=row[5]
        if (position%2) == 2:
            solventpos=str(solvents['MeOH2'])
        else:
            solventpos=str(solvents['MeOH'])
        pipettebig.aspirate(mvolume,solventsource[solventpos])
        pipettebig.dispense(mvolume,finplate.wells()[position])
        
        pipettebig.drop_tip()
        
    protocol.comment('End of Protocol.')