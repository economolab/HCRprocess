# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 15:52:09 2026

@author: jpv88
"""

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

# H1 contains I1* and I2
# H2 contains I1 and I2*

B1I1 = 'gAggAgggCAgCAAACgggAAgAgTCTTCCTTTACg'
B1I2 = 'gCATTCTTTCTTgAggAgggCAgCAAACgggAAgAg'
B2I1 = 'CCTCgTAAATCCTCATCAATCATCCAgTAAACCgCC'
B2I2 = 'AgCTCAgTCCATCCTCgTAAATCCTCATCAATCATC'
B3I1 = 'gTCCCTgCCTCTATATCTCCACTCAACTTTAACCCg'
B3I2 = 'AAAgTCTAATCCgTCCCTgCCTCTATATCTCCACTC'
B4I1 = 'CCTCAACCTACCTCCAACTCTCACCATATTCgCTTC'
B4I2 = 'CACATTTACAgACCTCAACCTACCTCCAACTCTCAC'
B5I1 = 'CTCACTCCCAATCTCTATCTACCCTACAAATCCAAT'
B5I2 = 'CACTTCATATCACTCACTCCCAATCTCTATCTACCC'

test = 'GAGCTCTGCCCTCATCAAAGATGGC'
testSeq = Seq(test)

ampSeqs = [B1I1, B1I2, B2I1, B2I2, B3I1, B3I2, B4I1, B4I2, B5I1, B5I2]
ampSeqs = [Seq(seq) for seq in ampSeqs]
ampSeqsTm = [mt.Tm_NN(seq, nn_table=mt.DNA_NN3) for seq in ampSeqs]

#%%

# 20X SSC is 3 M sodium chloride and 300 mM trisodium citrate so 3.9 M Na+
# 5X SSC is 0.975 M Na+ or 975 mM Na+
# 2X SSC is 0.39 M Na+ or 390 mM Na+
# PBS is 0.137 M Na+ and 0.0027 M K+

salt_correction_5X_SSC = mt.salt_correction(Na=975)
salt_correction_2X_SSC = mt.salt_correction(Na=390)
salt_correction_PBS = mt.salt_correction(Na=137, K=2.7)

fmd = 65
salt_correction_5X_SSC_fmd = mt.salt_correction(Na=975*((100-fmd)/100))
salt_correction_2X_SSC_fmd = mt.salt_correction(Na=390*((100-fmd)/100))
salt_correction_PBS_fmd = mt.salt_correction(Na=137*((100-fmd)/100), K=2.7*((100-fmd)/100))

#%%

import pandas as pd

df = pd.read_csv("probe_seqs.txt", 
                 sep="\t", 
                 header=None, 
                 names=['index', 'name', 'seq'])

probe_seqs = df['seq'].values
probe_seqs = [''.join(c for c in s if not c.islower()) for s in probe_seqs]

probeSeqsTm = [mt.Tm_NN(seq, 
                        nn_table=mt.DNA_NN3, 
                        Na=1000) for seq in probe_seqs]

probeSeqsTm = [mt.chem_correction(Tm, fmd=fmd) for Tm in probeSeqsTm]
df['Tm'] = probeSeqsTm

#%%

name = []
Tm = []

for i in range(0, len(df) - 1, 2):
    
    row1 = df.iloc[i]
    row2 = df.iloc[i+1]
    
    name.append("_".join(row1['name'].split("_")[:2]))
    Tm.append(min([row1['Tm'], row2['Tm']]))
    
df_min_Tm = pd.DataFrame({'name': name, 'Tm': Tm})
    
#%%

import pandas as pd

df = pd.read_csv("probe_seqs.txt", 
                 sep="\t", 
                 header=None, 
                 names=['index', 'name', 'seq'])

probe_seqs = df['seq'].values
probe_seqs = [''.join(c for c in s if not c.islower()) for s in probe_seqs]

probeSeqsTm = [mt.Tm_NN(seq, 
                      nn_table=mt.R_DNA_NN1, 
                      saltcorr=salt_correction_2X_SSC_fmd) for seq in probe_seqs]

newProbeSeqsTm = [mt.chem_correction(Tm,fmd=65) for Tm in probeSeqsTm]

df['Tm'] = newProbeSeqsTm

mask = ['Snap25' in s for s in df['name']]
