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

core_name= "DATA/Output_ScenA1_mini_time"
FA1   = np.array([0.01, 0.1])
CA1   = np.array([10,50,90])/100

EA1   = np.loadtxt(core_name + "E.txt", delimiter=',')
GA1   = np.loadtxt(core_name + "G.txt", delimiter=',')

arrSA1     = np.zeros((3, len(GA1), len(CA1), len(EA1), len(FA1)))
Req_St_A1  = np.zeros((len(GA1), len(CA1), len(EA1), len(FA1)))
Req_End_A1 = np.zeros((len(GA1), len(CA1), len(EA1), len(FA1)))
ci = 0
for c in CA1:
    fi = 0
    for f in FA1:
        name_file = core_name + "_f=" + "%03d"%((int(f*100)))+ "_c=" + "%03d"%((int(c*100)))
        arrSA1[0, :,ci , : , fi]= np.loadtxt(name_file + "_S0.txt", delimiter=',')
        arrSA1[1, :,ci , : , fi]= np.loadtxt(name_file + "_S1.txt", delimiter=',')
        arrSA1[2, :,ci , : , fi]= np.loadtxt(name_file + "_S2.txt", delimiter=',')
        
        Req_St_A1[:,ci , : , fi]= np.loadtxt(name_file + "_Time_G1_H0.txt", delimiter=',')
        Req_End_A1[:,ci , : , fi] = np.loadtxt(name_file + "_Time_G1_H1.txt", delimiter=',')
        
        fi += 1
    ci += 1
    
    
#%% get data A2
core_name= "DATA/Output_ScenA2_mini_time"
FA2   = np.array([0.01, 0.1])
CA2   = np.array([10**2, 10**4, 10**6])

EA2   = np.loadtxt(core_name + "E.txt", delimiter=',')
GA2   = np.loadtxt(core_name + "G.txt", delimiter=',')

arrSA2    = np.zeros((3, len(GA2), len(CA2), len(EA2), len(FA2)))
Req_St_A2  = np.zeros((len(GA2)  , len(CA2), len(EA2), len(FA2)))
Req_End_A2 = np.zeros((len(GA2)  , len(CA2), len(EA2), len(FA2)))

ci = 0
for c in CA2:
    fi = 0
    for f in FA2:
        name_file = core_name + "_f=" + "%03d"%((int(f*100)))+  "_c=" + "%03d"%(int(np.log10(c)))
        arrSA2[0, :,ci , : , fi]= np.loadtxt(name_file + "_S0.txt", delimiter=',')
        arrSA2[1, :,ci , : , fi]= np.loadtxt(name_file + "_S1.txt", delimiter=',')
        arrSA2[2, :,ci , : , fi]= np.loadtxt(name_file + "_S2.txt", delimiter=',')
        
        Req_St_A2[:,ci , : , fi]= np.loadtxt(name_file + "_Time_G1_H0.txt", delimiter=',')
        Req_End_A2[:,ci , : , fi] = np.loadtxt(name_file + "_Time_G1_H1.txt", delimiter=',')
        
        fi += 1
    ci += 1

#%% get data A2
# core_name= "DATA/Output_ScenB_mini_"
# FB   = np.array([0.01, 0.1])
# CB   = np.array([10**4, 10**6])

# EB   = np.loadtxt(core_name + "E.txt", delimiter=',')
# GB   = np.loadtxt(core_name + "G.txt", delimiter=',')

# arr_S_B    = np.zeros((3, len(GB), len(CB), len(EB), len(FB)))
# Req_St_B  = np.zeros((len(GB)  , len(CB), len(EB), len(FB)))
# Req_End_B = np.zeros((len(GB)  , len(CB), len(EB), len(FB)))

# ci = 0
# for c in CB:
#     fi = 0
#     for f in FB:
        
#         name_file = core_name + "_f=" + "%03d"%((int(f*100))) + "_c=" + "%03d"%(int(np.log10(c)))
#         arr_S_B[0, :,ci , : , fi]= np.loadtxt(name_file + "_S0.txt", delimiter=',')
#         arr_S_B[1, :,ci , : , fi]= np.loadtxt(name_file + "_S1.txt", delimiter=',')
#         arr_S_B[2, :,ci , : , fi]= np.loadtxt(name_file + "_S2.txt", delimiter=',')
        
#         Req_St_B= np.loadtxt(name_file + "_Time_G1_H0.txt", delimiter=',')
#         Req_End_B = np.loadtxt(name_file + "_Time_G1_H1.txt", delimiter=',')
        
#         fi += 1
#     ci += 1
#%% Cloud figures
markA1 ='o'
markA2 ='x'
markB ='.'
colorA1  =  "tab:green"#(0.5*0.7   , 0.4*0.7   , 1*0.7  ) # = 0.9 , 0.8
colorA2  =  "tab:blue"#(0         , 158/255   , 115/255) # = 2.1   , 0.6
colorB   =  "tab:orange"#(0.8       , 0.4       , 0.5    ) # = 1.9      , 0.9
legend_elements = [Line2D([0], [0], color=colorA1, alpha = 0.5, lw=4, label='A1'),
                   Line2D([0], [0], color=colorA2, alpha = 0.5, lw=4, label='A2'),
                   Line2D([0], [0], color=colorB , alpha = 0.5, lw=4, label='B')]

legend_elements2 = [Line2D([0], [0], color=colorA1, alpha = 0.5, lw=4, linestyle= ':', label='A1'),
                   Line2D([0], [0], color=colorA2, alpha = 0.5, lw=4, linestyle= '--',  label='A2'),
                   Line2D([0], [0], color=colorB , alpha = 0.5, lw=4, label='B')]

#%% PLots salinity

iif1=  0
iif2= -1

iic1=  -1
iic2= 0

fs  =  15

fig, ax = plt.subplots()
axins = ax.inset_axes([0.4, 0.47, 0.4, 0.4])
S= 1

for ie in [-1, 0]:
    y1 = arrSA1[S, :, iic1, ie, iif1] 
    y2 = arrSA1[S, :, iic2, ie, iif2]
    ax.fill_between(GA1, y1, y2, where= y2 <=350, interpolate=True, color='grey',
                 alpha=0.5)
    axins.fill_between(GA1, y1, y2, where= y2 <=350, interpolate=True, color='grey',
                 alpha=0.5)
    
    y1 = arrSA2[S, :, iic1, ie, iif1] 
    y2 = arrSA2[S, :, iic2, ie, iif2] 
    ax.fill_between(GA2, y1, y2, where= y2 <=350, interpolate=True, color='grey',
                 alpha=0.5)
    axins.fill_between(GA2, y1, y2, where= y2 <=350, interpolate=True, color='grey',
                 alpha=0.5)
    
    ax.plot(GA1, arrSA1[S, :, iic1, ie, iif1] , color = colorA1)
    ax.plot(GA1, arrSA1[S, :, iic2, ie, iif2] , color = colorA1)
    axins.plot(GA1, arrSA1[S, :, iif1, ie, iif1], color = colorA1)
    axins.plot(GA1, arrSA1[S, :, iic2, ie, iif2], color = colorA1)
    
    ax.plot(GA2, arrSA2[S, :, iic1, ie, iif1] , color = colorA2)
    ax.plot(GA2, arrSA2[S, :, iic2, ie, iif2] , color = colorA2)
    axins.plot(GA2, arrSA2[S, :, iic1, ie, iif1] , color = colorA2)
    axins.plot(GA2, arrSA2[S, :, iic2, ie, iif2] , color = colorA2)

    
    #ie+=1
ax.hlines(145   , GA1[0], GA1[-1],  linewidth= 2, linestyle = ':', color= 'k')
axins.hlines(145, GA1[0], GA1[-1],  linewidth= 2, linestyle = ':', color= 'k')
ax.hlines(350   , GA1[0], GA1[-1],  linewidth= 2, linestyle = '-', color= 'k')
axins.hlines(350, GA1[0], GA1[-1],  linewidth= 2, linestyle = '-', color= 'k')

# inset axes....

# sub region of the original image
x1, x2, y1, y2 =GA1[0], 0.6*10**3, 310, 350
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)


ax.indicate_inset_zoom(axins, edgecolor="black")
axins.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xscale('log')
#axins.set_xscale('log')


axins.set_xticklabels([])
#axins.set_yticklabels([ ])
axins.set_xticks([])
#axins.set_yticks([])
ax.set_xlim([GA1[0], GA1[-1]])
ax.set_ylim([35, 350])

ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_xlabel('low $\longleftarrow$ strait efficiency  $\longrightarrow$ high', color='black', fontsize= fs)

#plt.xlabel('Connection to Atlantic')
plt.ylabel('Salinity [$kg/m³$]', fontsize= fs)
plt.title("upper box", fontsize= fs*1.2)

plt.legend(handles=legend_elements, loc='upper right')
plt.tight_layout()
# fig.savefig( figname+ "UpperBox.png") 
# fig.savefig( figname+ "UpperBox.svg", format="svg") 
#plt.close() 


#%% PLOTS Time

iif1=  0
iif2= -1

iic1=  -1
iic2= 0

fs  =  15

fig, ax = plt.subplots()


for ie in [-1, 0]:
    y1 = Req_St_A1[:, iic1, ie, iif1]/1000
    y2 = Req_St_A1[:, iic2, ie, iif2]/1000
    ax.fill_between(GA1, y1, y2, where= y2 <=350, interpolate=True, color='grey',
                 alpha=0.5)

    
    y1 = Req_St_A2[:, iic1, ie, iif1]/1000
    y2 = Req_St_A2[:, iic2, ie, iif2]/1000
    ax.fill_between(GA2, y1, y2, where= y2 <=350, interpolate=True, color='grey',
                 alpha=0.5)

    
    ax.plot(    GA1, Req_St_A1[ :, iic1, ie, iif1]/1000 , color = colorA1)
    ax.plot(    GA1, Req_St_A1[ :, iic2, ie, iif2]/1000 , color = colorA1)

    
    ax.plot(    GA2, Req_St_A2[ :, iic1, ie, iif1]/1000 , color = colorA2)
    ax.plot(    GA2, Req_St_A2[ :, iic2, ie, iif2]/1000 , color = colorA2)

ax.set_xscale('log')
ax.set_xlim([GA1[0], GA1[-1]])


ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_xlabel('low $\longleftarrow$ strait efficiency  $\longrightarrow$ high', color='black', fontsize= fs)

#plt.xlabel('Connection to Atlantic')
plt.ylabel('time [kyr]', fontsize= fs)
plt.title("Time until first box reaches halite", fontsize= fs*1.2)

plt.legend(handles=legend_elements, loc='upper right')
plt.tight_layout()
# fig.savefig( figname+ "UpperBox.png") 
# fig.savefig( figname+ "UpperBox.svg", format="svg") 
#plt.close() 

#%% PLOTS Time

iif1=  0
iif2= -1

iic1=  -1
iic2= 0

fs  =  15

fig, ax = plt.subplots()


for ie in [-1, 0]:
    y1 = Req_End_A1[:, iic1, ie, iif1]/1000
    y2 = Req_End_A1[:, iic2, ie, iif2]/1000
    ax.fill_between(GA1, y1, y2, where= y2 <=350, interpolate=True, color='grey',
                 alpha=0.5)

    
    y1 = Req_End_A2[:, iic1, ie, iif1]/1000
    y2 = Req_End_A2[:, iic2, ie, iif2]/1000
    ax.fill_between(GA2, y1, y2, where= y2 <=350, interpolate=True, color='grey',
                 alpha=0.5)

    
    ax.plot(    GA1, Req_End_A1[ :, iic1, ie, iif1]/1000 , color = colorA1)
    ax.plot(    GA1, Req_End_A1[ :, iic2, ie, iif2]/1000 , color = colorA1)

    
    ax.plot(    GA2, Req_End_A2[ :, iic1, ie, iif1]/1000 , color = colorA2)
    ax.plot(    GA2, Req_End_A2[ :, iic2, ie, iif2]/1000 , color = colorA2)

ax.set_xscale('log')
ax.set_xlim([GA1[0], GA1[-1]])


ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_xlabel('low $\longleftarrow$ strait efficiency  $\longrightarrow$ high', color='black', fontsize= fs)

#plt.xlabel('Connection to Atlantic')
plt.ylabel('time [kyr]', fontsize= fs)
plt.title("Time until second box reaches halite", fontsize= fs*1.2)

plt.legend(handles=legend_elements, loc='upper right')
plt.tight_layout()
# fig.savefig( figname+ "UpperBox.png") 
# fig.savefig( figname+ "UpperBox.svg", format="svg") 
#plt.close() 