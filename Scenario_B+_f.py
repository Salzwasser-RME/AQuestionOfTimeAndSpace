# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 10:58:58 2022

@author: ronjaebner

Scenario B+
The Mediterranean Sea, as an invinite reservoir, is feeding a marginal basin.\
The MS is saturated in both gypsum and Halite
the Marginals Basin is experiencing evaporation and river input. 
The river do carry solved ions, but no chlorine

RQ: can we precipitate gypsum but no halite in this margin?

relevant equations


Problems:
for large g, dt has to be small otherwise the volume of the basin is too small
for the fluxes


chemical input Gaillardet et al. 1999, table 1
River= [NaCl, CaSo4] #kg/m3
Rhone= [0.03, 0.07]
Po   = [0.03, 0.09]
Nile = [0.07, 0.07]
Ebro = [0.12, 0.19]

f= (r-eA)/ fwb = fwb/R
f= [(c_r/ c_MS)h, (_r/ c_MS)g]
give area and evapo, calculate river inflow
R= (eA)/(1-f)

"""

#%% import necessary libraries and data
import           numpy   as np
import matplotlib.pyplot as plt
#import CMD_AP2_functions as fnc

    


#%% Define paramters
# time
T       = 20_000 #periode of cycle
t_max   = 1*T  # 
yr2sc   = 3600*24*365.25
dt      = 1  #[year]
i_max   = int(t_max/dt)
t_vec   = np.arange(0,i_max )

#Basin
A= 0.001*2.5*10**12
D=200
V= A*D #m³, random volume at that point
DTV= dt*yr2sc/V

g= 10**2
#%%
#%% Input rivers
#River= [NaCl, CaSo4] #kg/m3
River_str=['no ions',
            #'only CaSO4',
            'Rhone',
            'Nile',
            'Po',
            'Ebro']
Test0= [  0.00, 0.00,  0.00] # 0
Test1= [  0.00, 0.10,  0.00] # 1
Rhone= [  0.03, 0.07,  0.25] # 2
Po   = [  0.03, 0.09,  0.24] # 3
Nile = [  0.07, 0.07,  0.24] # 4
Ebro = [  0.12, 0.19,  0.21] # 5
MS   = [271.1 , 5.25, 72.90]

River_ions=np.array((Test0,Rhone, Nile, Po, Ebro))

e       = 0.5#m/yr
EP      = (A*e/yr2sc)
#factor  = np.logspace(-4.4*10**(-3), 0.017, 300)#
factor    = np.linspace(0.999, 1.04, 400)
R       = EP*factor

data= np.zeros([6, 5,len(R)])
        # MS
SM_H= MS[0]
SM_G= MS[1] #kg/m³
SM_R= MS[2]
SM  = 350
rr=-1
while rr <4:
    rr+=1
    RR=-1
    while RR < len(R)-1:
        RR+=1
        
        #River
        SR_H= River_ions[rr,0]
        SR_G= River_ions[rr,1] #kg/m³
        SR_R= River_ions[rr,2]
        SR = SR_H+ SR_G+ SR_R

        #% evolution over time inclusing
        SB = SM
        SB_G= SM_G
        SB_H= SM_H
        SB_R= SM_R
        fwb   = R[RR] - EP     #m3/s
        
        tt=-1
        while tt <t_max:
            tt+=1 
            # volume conservation    
            Q     = g * abs(SB - SM )   # m3/s
            F_in  = Q + max(0 ,(-fwb))  # m3/s
            F_out = Q + max(0 ,( fwb))  # m3/s
            
            # change in ions per dt
            SB_G_dt  = (SR_G*R[RR] + SM_G*F_in - SB_G*F_out)*DTV
            SB_H_dt  = (SR_H*R[RR] + SM_H*F_in - SB_H*F_out)*DTV
            SB_R_dt  = (SR_R*R[RR] + SM_R*F_in - SB_R*F_out)*DTV
            
            SB_G_tmp = SB_G + SB_G_dt
            SB_H_tmp = SB_H + SB_H_dt
            SB_R_tmp = SB_R + SB_R_dt
            
            #check for excess ions
            G_ex  = max(0, SB_G_tmp-SM_G)   
            H_ex  = max(0, SB_H_tmp-SM_H)
            
            # evolv salinity
            SB_G  = SB_G_tmp - G_ex
            SB_H  = SB_H_tmp - H_ex
            SB_R  = SB_R_tmp
            
            SB    = SB_G + SB_H + SB_R
            
        data[0,rr, RR] = SB_G
        data[1,rr, RR] = SB_H
        data[2,rr, RR] = SB_R
        
        data[3,rr, RR] = ((G_ex*V /dt)/2300)/A # kg/yr
        data[4,rr, RR] = ((H_ex*V /dt)/2200)/A # kg/yr
        data[5,rr, RR] = Q
#%% NICE pLOTS WITH THICKNESS OF HALITE ND GYPSUM
from matplotlib.patches import Patch
from matplotlib.legend import Legend
c_hal = "green"
c_gyp = "blue"
fig, ax = plt.subplots(5, 1, sharex='col', sharey='row', figsize=(9, 6), dpi=200)
ax[0].set_title('Gypsum and Halite precipitation for different rivers')
rr=-1
while rr <4:
    rr+=1 
    ax[rr].vlines(0,0, 0.05, color= "grey")
    ax[rr].scatter(0,0, color= "w", label= River_str[rr])
    ax[rr].fill_between(1-1/factor, data[3,rr,:]*1000,0, color= c_gyp, alpha = 0.5)
    ax[rr].fill_between(1-1/factor, data[4,rr,:]*1000,0, color= c_hal, alpha = 0.5)
    ax[rr].set_ylim([0, 0.05])
    ax[rr].legend(markerscale=1., scatterpoints=1, fontsize=10, loc = "center right")
    ax[rr].spines['left'].set_visible(False)
    ax[rr].spines['right'].set_visible(False)
    ax[rr].spines['top'].set_visible(False)
    
    #ax[rr].spines['bottom'].set_visible(False)
ax[rr].set_xlim([1-1/factor[0], 1-1/factor[-1]])

fig.text(0.04, 0.5, 'precipitation rate [m/kyr]', va='center', rotation='vertical')
legend_elements = [Patch(facecolor=c_gyp , edgecolor='w',label='Gypsum'),
                   Patch(facecolor=c_hal, edgecolor='w',label='Halite')]
leg = Legend(ax[0],legend_elements,["Gypsum", "Halite"],loc='lower center', frameon=False)
ax[0].add_artist(leg)
ax[4].set_xlabel( "fwb/R" )
plt.subplots_adjust(left=0.125,
                    bottom=0.125,
                    right=0.95,
                    top=0.9,
                    wspace=0.,
                    hspace=0.)
#%% NICE pLOTS WITH THICKNESS OF HALITE ND GYPSUM
from matplotlib.patches import Patch
from matplotlib.legend import Legend
c_hal = "green"
c_gyp = "blue"
fig, ax = plt.subplots(5, 1, sharex='col', sharey='row', figsize=(9, 6), dpi=200)
ax[0].set_title('Gypsum and Halite precipitation for different rivers',  fontsize= 16)
rr=-1
while rr <4:
    rr+=1 
    ax[rr].vlines(1,0, 0.05, color= "black", linestyle= ":" )
    ax[rr].scatter(0,0, color= "w", label= River_str[rr])
    ax[rr].fill_between(factor, data[3,rr,:]*1000,0, color= c_gyp, alpha = 0.5)
    ax[rr].fill_between(factor, data[4,rr,:]*1000,0, color= c_hal, alpha = 0.5)
    ax[rr].set_ylim([0.00, 0.05])
    ax[rr].set_xlim([0.99, 1.04])
    ax[rr].legend(markerscale=1., scatterpoints=1, fontsize=10, loc = "center right")
    ax[rr].spines['left'].set_visible(False)
    ax[rr].spines['right'].set_visible(False)
    ax[rr].spines['top'].set_visible(False)
    ax[rr].set_yticks([0, 0.04])
ax[rr].set_xlim([0.999, 1.04])#set_xlim([factor[0], factor[-1]])
ax[rr].set_xticks([1.0,  1.01 , 1.02, 1.03, 1.04])
fig.text(0.04, 0.5, 'precipitation after 1kyr [m]', va='center', rotation='vertical', fontsize= 14)
legend_elements = [Patch(facecolor=c_hal, edgecolor='w',label='Halite', alpha= 0.5),
                   Patch(facecolor=c_gyp , edgecolor='w',label='Gypsum', alpha= 0.5)]
leg = Legend(ax[0],legend_elements,["Halite", "Gypsum"],loc='lower center', frameon=False, fontsize= 14)
ax[0].add_artist(leg)
ax[4].set_xlabel( "R/EP" , fontsize= 14)
plt.subplots_adjust(left=0.125,
                    bottom=0.125,
                    right=0.95,
                    top=0.9,
                    wspace=0.,
                    hspace=0.1)
#%% Analysis
       # for every river ion composition we want to know if
        # data[3,rr, :]>0 && data[4,rr, :]==0
#%% linear Analysis
# This is a 1 dim analysis focussing on the  fwb
        # fwb= R-Ae
        # plot two lines for each river composition (halite prec, gypsum prec)
        # add calculated limits of R_test
        #color area where  data[3,rr, :]>0 && data[4,rr, :]==0 is true
# data_G= np.empty([len(R)]) 
# data_H= np.empty([len(R)]) 
       
# fig, axes = plt.subplots(6, 1, sharex='col', sharey='row', figsize=(9, 6), dpi=200)
# axes[0].set_title('Would Gypsum precipitate?')
# yticks=[0,1]
# rr=-1
# kw = dict(color='grey')
# while rr <5:
#     rr+=1
#     data_G= np.NaN
#     data_H= np.NaN
#     data_G= data[3,rr, :]
#     data_G[data_G>0]= 1
#     data_G[data_G!=1]= 0
#     data_H= data[4,rr, :]
#     data_H[data_H>0]=-1
#     data_H[data_H!=-1]=0

#     axes[rr].vlines(f[rr,0], 0,1)
#     axes[rr].vlines(f[rr,1], 0,1)
#     scatter= axes[rr].scatter(factor, data_H+data_G, label= River_str[rr],
#                                facecolors="black",
#                                marker= "+", s= 20)   
    
#     axes[rr].grid(color='w')
#     axes[rr].legend(markerscale=1., scatterpoints=1, fontsize=10)
#     plt.gray()
#     axes[rr].set_yticks(yticks)
#     axes[rr].set_yticklabels(['no', 'yes'])
#     axes[rr].set_ylim([-0.4, 1.4])
#     axes[rr].spines['left'].set_visible(False)
#     axes[rr].spines['right'].set_visible(False)
#     axes[rr].spines['top'].set_visible(False)
#     axes[rr].spines['bottom'].set_visible(False)
#     axes[rr].set_facecolor('lightgray' ) 
#     axes[rr].set_axisbelow(True)
# axes[5].set_xlabel('R/EP')
# dateiName= 'ScenB+_Comparing 6compositions_relativeFWB'
# plt.savefig(dateiName)
# plt.savefig(dateiName + '.eps', format='eps')

#    axes[rr-1].scatter(((A*e/yr2sc))/(R- (A*e/yr2sc) ), data_H+data_G)        
#    axes[rr-1].vlines((A*e/yr2sc)/(R_test[rr,0]-( A*e/yr2sc)), 0,1,'r')
#    axes[rr-1].vlines((A*e/yr2sc)/(R_test[rr,1]- (A*e/yr2sc)), 0,1)
        
#%% PLOTS OLD
#color_1='blue'
#color_2='green'
#color_3='orange'
#scriptSize= 16
#lineW=2
#indbeg=200 # first 200 timesteps are ignored
#f, axes = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(9, 6), dpi=200)
#distance = []        
#axes[0].plot(t_vec[2:-2]*dt , data[0,2:-2], linestyle='-', color= color_1, linewidth=lineW,
#    label='SG margin')
#axes[0].plot(t_vec[2:-2]*dt , data[1,2:-2], linestyle=':', color= color_2,linewidth=lineW,
#    label='SH margin')
#axes[0].plot(t_vec[2:-2]*dt , data[2,2:-2], linestyle=':', color= color_3,linewidth=lineW,
#    label='SR margin')
#axes[0].plot(t_vec[2:-2]*dt , data[0,2:-2]+data[1,2:-2]+data[2,2:-2], linestyle='-', color= 'k',linewidth=lineW,
#    label='SR margin')
#axes[0].grid()
#axes[0].set_ylabel('conc. [kg/m3]', fontsize=scriptSize)
#axes[0].legend(loc='best', fontsize=scriptSize-2) 
#axes[0].set_title('Basin fed by saturated MS,  fwb=0 \n g= ' 
#          + "{:.0E}".format(g) + ' Sv/(kg/m³),/n  mG  =' 
#          + str(int(data[3,-2]*1000))  + "g/kyr, mH = " + str(int(data[4,-2]*1000000))+ "mg/kyr")
#
#axes[1].plot(t_vec[2:-2]*dt,data[3,2:-2], color=color_1, linewidth=lineW,
#    label='Gypsum prec.')
#axes[1].plot(t_vec[2:-2]*dt,data[4,2:-2], color=color_2, linewidth=lineW,
#    label='Halite prec.')
#axes[1].grid()
#axes[1].set_ylabel('precip. [kg/yr]', fontsize=scriptSize)
#axes[1].set_xlabel('time [yr]', fontsize=scriptSize)
#axes[1].legend(loc='best', fontsize=scriptSize-2) 
##
##dateiName= scn + 'linesPlot'
##plt.savefig(dateiName)
##%% 
#test= data[0,:]+data[1,:]+data[2,:] -SM