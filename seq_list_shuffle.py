# RPA+FH model for single charge sequence
# Short-range interaction contributes to only the k=0 FH term
# The FH term ehs parameters follow the definition in the PRL paper

# ver Git.1 Apr 14, 2020
# Upload to github

# Jan 1, 2019
# - Add function for pH-dependent charges

import numpy as np
import pipit
from pipit import shuffle

R_pKa = 12.10
H_pKa = 6.04
K_pKa = 10.67
D_pKa = 3.71
E_pKa = 4.15

def qi_calc(x):
    return 10**x/(1+10**x)

def get_the_charge(seqname, pH=None):
    seq_noDE=polymers.EB1_linker[1:3]+polymers.EB1_linker[4]+polymers.EB1_linker[6:14]+polymers.EB1_linker[15:]
    shuffuled_noDE=shuffle.seq(seq_noDE, 1,66)
    shuffuled_addDE='D'+shuffuled_noDE[0:2]+'D'+shuffuled_noDE[2:3]+'D'+shuffuled_noDE[3:11]+'E'+shuffuled_noDE[11:]
    
    #shuffuled_linker=shuffle.seq(polymers.EB1_linker, 1,71)
    #the_seq = shuffuled_linker+polymers.EB1_coil+polymers.EB1_tail
    the_seq = shuffuled_addDE+polymers.EB1_coil+polymers.EB1_tail
    N = len(the_seq)
    sigmai = np.zeros(N)
    use_pKa = False if pH == None else True
     
    for i in range(0, N):
        if the_seq[i] == 'D' :
            sigmai[i] = -qi_calc(pH-D_pKa) if use_pKa else -1
        elif the_seq[i] == 'E' :
            sigmai[i] = -qi_calc(pH-E_pKa) if use_pKa else -1
        elif the_seq[i] == 'R' :
            sigmai[i] = qi_calc(R_pKa-pH) if use_pKa else 1
        elif the_seq[i] == 'K' :
            sigmai[i] = qi_calc(K_pKa-pH) if use_pKa else 1
        elif the_seq[i] == 'H' :
            sigmai[i] = qi_calc(H_pKa-pH) if use_pKa else 0 
            if 'pH1' in seqname:
                print('seq is pH1')
                sigmai[i] = 1     
        else:
            sigmai[i] = 0

    return sigmai, N, the_seq




class polymers:

    EB1_coil = "QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ" 

    EB1_tail = "TDEGFVIPDE" + "GGPQEEQEEY"
    
    EB1_linker = \
    "D" + "GKDYDPVAAR" + "QGQETAVAPS" + "LVAPALNKPK" + "KPLTSSSAAP" +\
    "QRPISTQRTA" + "AAPKAGPGVV" + "RKNPGVGNG"


