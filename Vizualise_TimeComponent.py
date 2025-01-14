# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 14:28:50 2024

@author: rme_w
"""

import           numpy   as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker
#%%

'''
Load the data separately, this is not yet implemented
'''
#%% get data A1
SH = 350
SG = 145

run = "_mini_time_focus"
figname = "../qoTaS_FIGURE/Output" + run
#%%
core_name= "../qoTaS_DATA/Output_ScenA1" + run
FA1   = np.array([0.01, 0.1])
CA1   = np.array([10,50,90])/100

E   = np.loadtxt(core_name + "E.txt", delimiter=',')
G   = np.loadtxt(core_name + "G.txt", delimiter=',')

arrS      = np.zeros((3, len(G), len(CA1), len(E), len(FA1), 3))
T_First_G = np.zeros((len(G), len(CA1), len(E), len(FA1), 3))
T_First_H = np.zeros((len(G), len(CA1), len(E), len(FA1), 3))
T_All_G= np.zeros((len(G), len(CA1), len(E), len(FA1), 3))
T_All_H= np.zeros((len(G), len(CA1), len(E), len(FA1), 3))

conf = 0
ci = 0
for c in CA1:
    fi = 0
    for f in FA1:
        name_file = core_name + "_f=" + "%03d"%((int(f*100)))+ "_c=" + "%03d"%((int(c*100)))
        arrS[0, :,ci , : , fi, conf]= np.loadtxt(name_file + "_S0.txt", delimiter=',')
        arrS[1, :,ci , : , fi, conf]= np.loadtxt(name_file + "_S1.txt", delimiter=',')
        arrS[2, :,ci , : , fi, conf]= np.loadtxt(name_file + "_S2.txt", delimiter=',')
        
        T_First_G[:,ci , : , fi, conf] = np.loadtxt(name_file + "_time_first_gypsum.txt", delimiter=',')
        T_First_H[:,ci , : , fi, conf] = np.loadtxt(name_file + "_time_first_halite.txt", delimiter=',')
        T_All_G[:,ci , : , fi, conf]   = np.loadtxt(name_file + "_time_all_gypsum.txt"  , delimiter=',')
        T_All_H[:,ci , : , fi, conf]   = np.loadtxt(name_file + "_time_all_halite.txt"  , delimiter=',')
        
        fi += 1
    ci += 1
    
    
#%% get data A2
core_name= "../qoTaS_DATA/Output_ScenA2" + run
FA2   = np.array([0.01, 0.1])
CA2   = np.array([10**2, 10**4, 10**6])

conf = 1
ci = 0
for c in CA2:
    fi = 0
    for f in FA2:
        name_file = core_name + "_f=" + "%03d"%((int(f*100)))+  "_c=" + "%03d"%(int(np.log10(c)))
        arrS[0, :,ci , : , fi, conf]= np.loadtxt(name_file + "_S0.txt", delimiter=',')
        arrS[1, :,ci , : , fi, conf]= np.loadtxt(name_file + "_S1.txt", delimiter=',')
        arrS[2, :,ci , : , fi, conf]= np.loadtxt(name_file + "_S2.txt", delimiter=',')
        
        T_First_G[:,ci , : , fi, conf] = np.loadtxt(name_file + "_time_first_gypsum.txt", delimiter=',')
        T_First_H[:,ci , : , fi, conf] = np.loadtxt(name_file + "_time_first_halite.txt", delimiter=',')
        T_All_G[:,ci , : , fi, conf]   = np.loadtxt(name_file + "_time_all_gypsum.txt"  , delimiter=',')
        T_All_H[:,ci , : , fi, conf]   = np.loadtxt(name_file + "_time_all_halite.txt"  , delimiter=',')
        
        fi += 1
    ci += 1

#%% get data A2
core_name= "../qoTaS_DATA/Output_ScenB" + run
FB   = np.array([0.01, 0.1])
# CB   = np.array([10**2, 10**4, 10**6])
CB   = np.array([10**4, 10**6])
conf = 2
ci = 0
for c in CB:
    fi = 0
    for f in FB:
        
        name_file = core_name + "_f=" + "%03d"%((int(f*100))) + "_c=" + "%03d"%(int(np.log10(c)))
        arrS[0, :,ci , : , fi, conf]= np.loadtxt(name_file + "_S0.txt", delimiter=',')
        arrS[1, :,ci , : , fi, conf]= np.loadtxt(name_file + "_S1.txt", delimiter=',')
        arrS[2, :,ci , : , fi, conf]= np.loadtxt(name_file + "_S2.txt", delimiter=',')
        
        T_First_G[:,ci , : , fi, conf] = np.loadtxt(name_file + "_time_first_gypsum.txt", delimiter=',')
        T_First_H[:,ci , : , fi, conf] = np.loadtxt(name_file + "_time_first_halite.txt", delimiter=',')
        T_All_G[:,ci , : , fi, conf]   = np.loadtxt(name_file + "_time_all_gypsum.txt"  , delimiter=',')
        T_All_H[:,ci , : , fi, conf]   = np.loadtxt(name_file + "_time_all_halite.txt"  , delimiter=',')
        
        
        fi += 1
    ci += 1
#%% get rid of 0s
arrS[       arrS        ==0.0] = np.nan    
T_First_G[  T_First_G   ==0.0] = np.nan 
T_First_H[  T_First_H   ==0.0] = np.nan 
T_All_G[    T_All_G     ==0.0] = np.nan 
T_All_H[    T_All_H     ==0.0] = np.nan 
#%% Definitions for plots
markA1 ='o'
markA2 ='x'
markB ='.'
colorA1  =  "tab:green"
colorA2  =  "tab:blue"
colorB   =  "tab:orange"

legend_elements_3 = [Line2D([0], [0], color='k', alpha = 0.5, lw=3,linestyle= '-', label='Gypsum'),
                     Line2D([0], [0], color='k', alpha = 0.5, lw=3,linestyle= ':', label='Halite'),
                     Line2D([0], [0], color='w', alpha = 0.0, lw=1,linestyle= ':', label=' '),
                     Line2D([0], [0], color='k' , alpha = 0.5, lw=2              , label='e =0.25 m/yr'),
                     Line2D([0], [0], color='k' , alpha = 0.5, lw=3              , label='e =0.75 m/yr'),
                     Line2D([0], [0], color='w', alpha = 0.0, lw=1,linestyle= ':', label=' '),
                     Line2D([0], [0], color=colorA1, alpha = 0.5, lw=3, label='A1'),
                     Line2D([0], [0], color=colorA2, alpha = 0.5, lw=3, label='A2'),
                     Line2D([0], [0], color=colorB , alpha = 0.5, lw=3, label='B')]

legend_elements_4 = [Line2D([0], [0], color='k' , alpha = 0.5, lw=2              , label='e =0.25 m/yr'),
                     Line2D([0], [0], color='k' , alpha = 0.5, lw=3              , label='e =0.75 m/yr'),
                     Line2D([0], [0], color='w', alpha = 0.0, lw=1,linestyle= ':', label=' '),
                     Line2D([0], [0], color=colorA1, alpha = 0.5, lw=3, label='A1'),
                     Line2D([0], [0], color=colorA2, alpha = 0.5, lw=3, label='A2'),
                     Line2D([0], [0], color=colorB , alpha = 0.5, lw=3, label='B')]

color_conf = [colorA1, colorA2, colorB]

y_min = 0
y_max = 130


def get_last_g(conf, iie, iic, T_H, T_G, case):
    if case == 'first':
        y1 = T_G[ :, iic1, iie, iif1, conf]/1000
    elif case == 'met':
        y1 = ((T_All_H[ :, iic, iie, iif1, conf] - 
               np.maximum(T_First_H[ :, iic, iie, iif1, conf], 
                          T_All_G[ :, iic, iie, iif1, conf])   )/1000)
    else:
        print('Case must be either "first" or "met"')
        
    y1[y1==0.0]= np.nan
    iig = np.where(np.isnan(y1))
    return iig
 
#%% PLot time ti first gypsum and first halite for different configurations
iif1=  0
iif2= -1

iic1=  1
iic2= 1

fs  =  15

fig, ax = plt.subplots(dpi = 300)
S= 1
for conf in [0,1,2]:
    for iie in [1, 5]:
        if iie == 1:
            lw = 1
        elif iie ==5:
            lw = 2
        y1 = T_First_G[ :, iic1, iie, iif1, conf] /1000
        y2 = T_First_G[ :, iic2, iie, iif2, conf] /1000

        ax.plot(G, y1 , color = color_conf[conf], linestyle = '-', linewidth = lw)
        ax.plot(G, y2 , color = color_conf[conf], linestyle = '-', linewidth = lw)
        
        y1 = T_First_H[ :, iic1, iie, iif1, conf] /1000
        y2 = T_First_H[ :, iic2, iie, iif2, conf] /1000
        
        ax.plot(G, y1 , color = color_conf[conf], linestyle = ':', linewidth = lw)
        ax.plot(G, y2 , color = color_conf[conf], linestyle = ':', linewidth = lw)


conf = 1
iie = 5
iic = 1
iig = get_last_g(conf, iie, iic, T_First_H, T_First_G, 'first')
ax.vlines(G[iig[0][0]-1], y_min, y_max , color = 'gray', alpha = 0.5)


ax.set_xscale('log')
ax.set_xlim([10**2, 10**6])
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_xlabel('low $\longleftarrow$ strait efficiency  $\longrightarrow$ high', color='black', fontsize= fs)

plt.ylabel('time [$kyr$]', fontsize= fs)
plt.title("time to first gypsum and first halite", fontsize= fs*1.2)
plt.ylim([y_min, y_max])
plt.legend(handles=legend_elements_3, loc='upper right')
plt.tight_layout()
fig.savefig( figname+ "time_first_both.png") 
fig.savefig( figname+ "time_first_both.svg", format="svg") 


#%% PLot time the conditions take

iif1=  0
iif2= -1

iic1=  0
iic2=  -1

fs  =  15

fig, ax = plt.subplots(dpi = 300)
axins = ax.inset_axes([0.25, 0.1, 0.35, 0.8])
S= 1
for conf in [0,1,2]:
    for iie in [1, 5]:
        for iic in [0, -1]:
            if iie == 1:
                lw = 1
            elif iie ==5:
                lw = 2
            y1 = ((T_All_H[ :, iic, iie, iif1, conf] - 
                   np.maximum(T_First_H[ :, iic, iie, iif1, conf], 
                              T_All_G[ :, iic, iie, iif1, conf])   )/1000)
            
            y1[y1==0.0]= np.nan
            
            ax.plot(G, y1 , color = color_conf[conf], linewidth = lw)
            axins.plot(G, y1 , color = color_conf[conf], linewidth = lw)
            

conf = 0
iie = -1
iic = 1
iig = get_last_g(conf, iie, iic, T_All_H, T_All_G, 'met')
ax.vlines(G[iig[0][0]-1], y_min, y_max , color = 'gray', alpha = 0.5)
axins.vlines(G[iig[0][0]-1], y_min, y_max , color = 'gray', alpha = 0.5)
axins.plot(G, y1 , color = color_conf[conf], alpha = 0.5)

ax.set_xscale('log')
ax.set_xlim([10**2, 10**6])
x1, x2, y1, y2 =G[0], 0.4*10**3, 10**(-2), 50
axins.set_yscale('log')
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
ax.indicate_inset_zoom(axins, edgecolor="black")
axins.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_xlabel('low $\longleftarrow$ strait efficiency  $\longrightarrow$ high', color='black', fontsize= fs)

plt.ylabel('time [$kyr$]', fontsize= fs)
plt.title("duration of simultaneous precipitation", fontsize= fs*1.2)

plt.legend(handles=legend_elements_4, loc='upper right')
plt.ylim([y_min, y_max])
plt.tight_layout()
fig.savefig( figname+ "duration_conditions_met_extra.png") 
fig.savefig( figname+ "duration_conditions_met_extra.svg", format="svg") 
