#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:21:19 2022

@author: Ronja Ebner

Script in support of the paper 
_______________________________________________________________________________
A question of time and space: 
    A model approach to the synchronicity of gypsum and halite during the 
    Messinian Salinity Crisis
Ebner&Meijer, 2024
_______________________________________________________________________________

Configuration A1

This script calculates the salinties for the different boxes for the A1 
configuration for multiple values of e,g,f,c.


It's output describes the timespan it takes for the marginal box to reach 
halite saturation after the open box has reached gypsum
"""

import           numpy   as np

#%% time and stuff
print( "initiate" )
yr2sc = 24*60*60*365.25
dt    = 0.5


#%% 

F   = np.array([0.01, 0.1])
A0  = 2.5*10**12
D0  = 1500
Dint=  500

SG  = 145
SA  =  36
SH  = 350
#
G   =np.logspace(2, 4, num=900)
#
E_step = 12.5
Etmp   = np.arange(12.5, 80, E_step)
E      = Etmp/(100*yr2sc)
#scenario A1 (convection)
Ctmp   = np.array([10,50,90])
C      = Ctmp/(100)

kappa_mix   = 1*10**(-4)
kappa_conv  = 1*10**(-1)
d_mix       = 0.5*D0

arrSA1 = np.zeros((3, len(G), len(C), len(E), len(F)))
arrFA1 = np.zeros((2, len(G), len(C), len(E), len(F)))
arrTime = np.zeros((4, len(G), len(C), len(E), len(F)))
#%% loop
print( "start loop" )
fi=0
for f in F:
    V    = np.array([f*A0*Dint,(1-f)*A0*Dint, A0*(D0-Dint) ])
    ic=0
    for c in C:
        ie =0
        for e in E:
            ig = 0
            for g in  G:
                # t-max increases with decreasing g to
                # 1. make sure that equlibrium is reached (small g)
                # 2. minimise the runtime (large g)
                t_max = 1_000_000/(np.log10(g))**2
                S    = np.ones((3, int(t_max/dt)+1))*SA
                t=0
                ii=0
                while t< t_max-dt:
                    t+=dt
                    
                    Q= g*np.sqrt(S[1,ii]-SA)
                    Qin= Q+ e*A0
                    
                    F02 = (S[0, ii]>S[2, ii])*kappa_conv*f*A0*(S[0,ii]-S[2,ii])/(S[2,ii]*d_mix)
                    F20 = (1-c)*F02
                    F21 = c*F02
                    F10 = c*F02 + e*A0*f
                    mix = kappa_mix*A0*(1-f)*(S[2, ii]-S[1, ii])/d_mix
            
                    S[0,ii+1]= min(SH, S[0,ii] + (F10*S[1,ii] + F20*S[2,ii] 
                                          - F02*S[0,ii])*dt*yr2sc/V[0])
                    S[1,ii+1]= min(SH, S[1,ii] + (Qin*SA + F21*S[2, ii] 
                                          -(F10 + Q)*S[1, ii]  
                                          + mix)*dt*yr2sc/V[1])
                    S[2,ii+1]= min(SH, S[2,ii] + (F02*S[0,ii] 
                                          - (F21 + F20)*S[2,ii] 
                                          - mix)*dt*yr2sc/V[2])
                    ii+=1
                time_first_gypsum   = 0    
                time_all_gypsum     = 0
                time_first_halite   = 0
                time_all_halite     = 0
                        
                if max(S[[0,1],ii-1])>= SG: # at least one surface box has reached Gypsum  
                    index               = np.where(np.any(S[[0,1],:]>=SG, axis=0))                
                    time_first_gypsum   = index[0][0]*dt
                    index01             = np.where(S[[0,1],:]>=SG)
                    test                = np.rec.find_duplicate(index01[1].tolist())
                    if len(test)>0 :# check if both boxes have reached gypsum               
                        time_all_gypsum = test[0]*dt
                        
                if max(S[[0,1],ii-1])>= SH: # at least one surface box has reached Gypsum  
                    index               = np.where(np.any(S[[0,1],:]>=SH, axis=0))                
                    time_first_halite   = index[0][0]*dt
                    index01             = np.where(S[[0,1],:]>=SH)
                    test                = np.rec.find_duplicate(index01[1].tolist())
                    if len(test)>0 :# check if both boxes have reached gypsum               
                        time_all_halite = test[0]*dt

                arrTime[:, ig, ic, ie, fi] = [time_first_gypsum,
                                              time_all_gypsum ,
                                              time_first_halite,
                                              time_all_halite]
                arrSA1[:, ig, ic, ie, fi]= [S[0,ii], S[1,ii], S[2,ii]]
                arrFA1[:, ig, ic, ie, fi]= [Q,F02]
                
                ig+=1
            ie+=1
        ic+=1
    fi+=1

            
#%% Save arrays in txt files
print( "begin saving" )
name_dir="../qoTaS_DATA/Output_ScenA1_mini_time_focus"
ci=0
for c in C:
    fi=0
    for f in F:
        #  %5.4f' %(3.141592))
        name_file = name_dir+"_f=" + "%03d"%((int(f*100))) + "_c=" + "%03d"%(int(c*100)) #%(int(np.log10(c)))) 
        
        #S0
        data = np.squeeze(arrSA1[0,:,ci, :, fi])
        np.savetxt(name_file+"_S0.txt",data,delimiter=",")
        
        #S1
        data = np.squeeze(arrSA1[1,:,ci, :, fi])
        np.savetxt(name_file+"_S1.txt",data,delimiter=",")
        
        #S2
        data = np.squeeze(arrSA1[2,:,ci, :, fi])
        np.savetxt(name_file+"_S2.txt",data,delimiter=",")
        
        #Q
        data = np.squeeze(arrFA1[0,:,ci, :, fi])
        np.savetxt(name_file+"_Q.txt",data,delimiter=",")
        
        #Q
        data = np.squeeze(arrFA1[1,:,ci, :, fi])
        np.savetxt(name_file+"_F.txt",data,delimiter=",")
        
        #time_first_gypsum 
        data = np.squeeze(arrTime[0,:,ci, :, fi])
        np.savetxt(name_file+"_time_first_gypsum.txt",data,delimiter=",")
        
        #time_all_gypsum 
        data = np.squeeze(arrTime[1,:,ci, :, fi])
        np.savetxt(name_file+"_time_all_gypsum.txt",data,delimiter=",")
        
        #time_first_halite
        data = np.squeeze(arrTime[2,:,ci, :, fi])
        np.savetxt(name_file+"_time_first_halite.txt",data,delimiter=",")
        
        #time_all_halite
        data = np.squeeze(arrTime[3,:,ci, :, fi])
        np.savetxt(name_file+"_time_all_halite.txt",data,delimiter=",")
        fi+=1
    ci+=1

#E
data = np.squeeze(Etmp)
np.savetxt(name_dir+"E.txt",data,delimiter=",")
#G
data = np.squeeze(G)
np.savetxt(name_dir+"G.txt",data,delimiter=",")
print( "saving successful" )
