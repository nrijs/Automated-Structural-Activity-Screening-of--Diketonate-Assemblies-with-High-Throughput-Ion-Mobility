# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 12:23:36 2021
@author: oscar
"""

## The dilution step after the mixing step for the worklow found in Automated Structural Activity Screening of β-Diketonate Assemblies with High-Throughput Ion Mobility - Mass Spectrometry.

from opentrons import protocol_api
metadata = {'apiLevel':'2.7',
           'author':'Oscar Lloyd Williams'}

import csv
import subprocess

##pasting in CSV again
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
##reading csv data again
csv_data = csv_copy.splitlines()[1:] # Discard the blank first line.
csv_reader = csv.reader(csv_data)
csv_list=list(csv_reader)
  
           ##opentrons protocol
def run(protocol: protocol_api.ProtocolContext):
    
    
    
    
    # Defining the methanol solvent location
    solventsource=protocol.load_labware('olw_2schott_100000ul',8)
    solvents={'MeOH':'A1', 'MeOH2':'A2'}
    
    #Defining the locations of the HPLC vial containing plates for the final solutions
    finplate1 = protocol.load_labware('olw_40hplc_2000ul', 2)
    finplate2 = protocol.load_labware('olw_40hplc_2000ul', 3)
    
    #Defining the locations of the stock vials
    startplate1 = protocol.load_labware('olw_40hplc_2000ul', 5)
    startplate2 = protocol.load_labware('olw_40hplc_2000ul', 6)
    
    #Defining the location of the tipracks for two different pipettes
    tiprackbig = protocol.load_labware('opentrons_96_tiprack_1000ul', 4)
    tipracklil= protocol.load_labware('opentrons_96_tiprack_300ul', 9)
    
    #Defining the location of the pH testing location.
    phplate= protocol.load_labware('corning_96_wellplate_360ul_flat',1)
    
    
    #Defining the pipettes used and adjusting the clearance to minimise chance of error with incosistent prints.
    pipettebig = protocol.load_instrument('p1000_single_gen2', 'left',tip_racks=[tiprackbig])
    pipettesmol = protocol.load_instrument('p300_single_gen2', 'right',tip_racks=[tipracklil])
    pipettebig.well_bottom_clearance.dispense = 3
    
    pipettebig.pick_up_tip()
    

## depositing methanol into evey vial
    for row in csv_list:
        finplate=finplate1
        position= int(row[0])
        if position >39:
            position=position-40
            finplate=finplate2
        mvolume=950
        if (position%2) == 2:
            solventpos=str(solvents['MeOH2'])
        else:
            solventpos=str(solvents['MeOH'])
        pipettebig.aspirate(mvolume,solventsource[solventpos])
        pipettebig.dispense(mvolume,finplate.wells()[position])
        
    pipettebig.drop_tip()
     ## adding 50 uL of sample solution to each vial and adding 20uL to the pH plate     
    for row in csv_list:
        finplate=finplate1
        pipettesmol.pick_up_tip()
        positionog= int(row[0])
        position=positionog
        if positionog >39:
            position=positionog-40
            finplate=finplate2
        startplate=startplate1
        positionog= int(row[0])
        position=positionog
        if positionog >39:
            position=positionog-40
            startplate=startplate2
        pipettesmol.aspirate(70,startplate.wells()[position])
        pipettesmol.dispense(50,finplate.wells()[position])
        pipettesmol.dispense(20,phplate.wells()[positionog])
        pipettesmol.drop_tip()
    
    protocol.comment('End of Protocol.')