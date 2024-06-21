#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:21:19 2022

@author: ron

C-Question

Scenario A2

# for ploting. it is important to find the right balance between 
"""

import           numpy   as np

#%% time and stuff
yr2sc = 24*60*60*365.25

dt = 0.5
DT =dt*yr2sc
t_max = 40_000

#%% 

F   = np.array([0.01, 0.1])
A0  = 2.5*10**12
D0  = 1500
Dint=  500

SA  =  36
SH  = 350
#
G   =np.logspace(2, 3, num=30)
#
E_step = 0.05
Etmp   = np.arange(50, 110, E_step)
E      = Etmp/(100*yr2sc)
#
#scenario A2 (restricted margin)
C      = np.array([10**2, 10**4, 10**6])

kappa_mix   = 1*10**(-4)
kappa_conv  = 1*10**(-1)
d_mix       = 0.5*D0




S    = np.ones((3, int(t_max/dt)))*SA

arrSA2 = np.zeros((3, len(G), len(C), len(E), len(F)))
arrFA2 = np.zeros((2, len(G), len(C), len(E), len(F)))

#%% loop
fi=0
for f in F:
    V    = np.array([f*A0*Dint,(1-f)*A0*Dint, A0*(D0-Dint) ])
    ic=0
    for A2c in C:
        ie =0
        for e in E:
            ig = 0
            for g in  G:
                #%%
                # Scenario A2
                t_max = 1_000_000/(np.log10(g))**2
                S    = np.ones((3, int(t_max/dt)+1))*SA
                t=0
                ii=0
                while t< t_max-dt:
                    # Scenario A2
                    t+=dt
            
                    Q= g*np.sqrt(S[1,ii]-SA)
                    Qin= Q+ e*A0
            
                    F02 = (S[0,ii]>S[1,ii])*(S[0,ii]>S[2,ii])*A2c*(S[0,ii]  - S[2,ii])
                    F01 = (S[0,ii]>S[1,ii])*(S[0,ii]<S[2,ii])*A2c*(S[0,ii]  - S[1,ii])
                    F10 = F02 + F01 +e*A0*f
                    F21 = F02
                    
                    mix = kappa_mix*A0*(1-f)*(S[2, ii]-S[1, ii])/d_mix
            
                    S[0, ii+1]= min(SH, S[0,ii] 
                                        + (F10*S[1,ii] 
                                           - (F01 +F02)*S[0,ii])
                                                *DT/V[0])
                    S[1, ii+1]= min(SH, S[1,ii] 
                                        + (Qin*SA + F21*S[2,ii] + F01*S[0,ii] 
                                           - (F10+Q)*S[1, ii] + mix)
                                                *DT/V[1])
                    S[2, ii+1]= min(SH, S[2,ii] 
                                        + (F02*S[0,ii] 
                                           - F21*S[2,ii] -mix)
                                                *DT/V[2])
                    ii+=1
                    
                arrSA2[:, ig, ic, ie, fi]= [S[0,ii], S[1,ii], S[2,ii]]
                arrFA2[:, ig, ic, ie, fi]= [Q,F02]
                
                ig+=1
            ie+=1
        ic+=1
    fi+=1
            
#%% Save arrays in txt files
name_dir="DATA/Output_ScenA2_mini_"
ci=0
for c in C:
    fi=0
    for f in F:
        #  %5.4f' %(3.141592))
        name_file = name_dir+"_f=" + "%03d"%((int(f*100))) + "_c=" + "%03d"%(int(np.log10(c))) #%(int(np.log10(c)))) 
        
        #S0
        data = np.squeeze(arrSA2[0,:,ci, :, fi])
        np.savetxt(name_file+"_S0.txt",data,delimiter=",")
        
        #S1
        data = np.squeeze(arrSA2[1,:,ci, :, fi])
        np.savetxt(name_file+"_S1.txt",data,delimiter=",")
        
        #S2
        data = np.squeeze(arrSA2[2,:,ci, :, fi])
        np.savetxt(name_file+"_S2.txt",data,delimiter=",")
        
        #Q
        data = np.squeeze(arrFA2[0,:,ci, :, fi])
        np.savetxt(name_file+"_Q.txt",data,delimiter=",")
        
        #Q
        data = np.squeeze(arrFA2[1,:,ci, :, fi])
        np.savetxt(name_file+"_F.txt",data,delimiter=",")
        
        fi+=1
    ci+=1

#E
data = np.squeeze(Etmp)
np.savetxt(name_dir+"E.txt",data,delimiter=",")
#G
data = np.squeeze(G)
np.savetxt(name_dir+"G.txt",data,delimiter=",")
