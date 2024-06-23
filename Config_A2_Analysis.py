# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 11:49:00 2023

@author: rme_w

curve fitting the output to get rid of the jumps

#example
1. define function to fit on
def model_f(x,a,b,c):
  return a*(x-b)**2+c

# a, b and c as p0=[3,2,-16].

2. call fitting function
and feed with base function, data and parameters
popt, pcov = curve_fit(model_f, x_data, y_data, p0=[3,2,-16])

3. get the optimized parameters

4. plot the comparinson

a_opt, b_opt, c_opt = popt
x_model = np.linspace(min(x_data), max(y_data), 100)
y_model = model_f(x_model, a_opt, b_opt, c_opt) 
 
plt.scatter(x_data, y_data)
plt.plot(x_model, y_model, color='r')
plt.show()


"""
import           numpy   as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit

#%% Read the data from files
# Output_ScanA1_f=xxx_c=xxx_VAR.txt
core_name= "Data/minirun/DATA/Output_ScenA2_"
VAR = ['S0', 'S1', 'S2', 'F', 'Q']
F   = np.array([0.01, 0.1])
C   = np.array([10**2, 10**4, 10**6])
# Get the axis
E   = np.loadtxt(core_name + "E.txt", delimiter=',')
G   = np.loadtxt(core_name + "G.txt", delimiter=',')
# define arrays
S0  = np.zeros((len(G), len(C), len(E), len(F)))
S1  = np.zeros((len(G), len(C), len(E), len(F)))
S2  = np.zeros((len(G), len(C), len(E), len(F)))
FLUX= np.zeros((len(G), len(C), len(E), len(F)))
Q   = np.zeros((len(G), len(C), len(E), len(F)))

# get ze data from ze filez
ci = 0
for c in C:
    fi = 0
    for f in F:
        name_file = core_name + "_f=" + "%03d"%((int(f*100)))+ "_c=" + "%03d"%((int(np.log10(c))))
        S0[:, ci , : , fi]= np.loadtxt(name_file + "_S0.txt", delimiter=',')
        S1[:, ci , : , fi]= np.loadtxt(name_file + "_S1.txt", delimiter=',')
        S2[:, ci , : , fi]= np.loadtxt(name_file + "_S2.txt", delimiter=',')
        FLUX[:, ci , : , fi]= np.loadtxt(name_file + "_F.txt", delimiter=',')
        Q[:, ci , : , fi]= np.loadtxt(name_file + "_Q.txt", delimiter=',')
        fi += 1
    ci += 1
fi =0
ci= 0   
#%%
test= np.zeros((len(G), len(E)))
name_file = core_name + "_f=" + "%03d"%((int(f*100)))+ "_c=" + "%03d"%((int(np.log10(c))))
test=np.loadtxt(name_file + "_F.txt", delimiter=',')
#%% 
A0  = 2.5*10**12
D0  = 1500
Dint=  500
yr2sc=60*60*24*365.25
SA  =  36
SH  = 350
#%% Analysis

# Find Conditions when coeval precipitation would occur
# S0=Sh & S1<SH
LowerLim =  np.zeros((len(G), len(C), len(E), len(F)))
UpperLim =  np.zeros((len(G), len(C), len(E), len(F)))
arr_dS   =  np.zeros((len(G), len(C), len(E), len(F)))

LowerLim[:,:,:] = S0[:, :, :, :]
UpperLim[:,:,:] = S1[:, :, :, :]

UpperLim[UpperLim==SH]=np.nan
UpperLim[LowerLim<SH] =np.nan

LowerLim[LowerLim<SH]=np.nan
LowerLim[LowerLim==SH]=1

arr_dS[:,:,:,:]= SH*LowerLim[:,:,:,:]-UpperLim[:,:,:, :]
#%% Analysis
LowerLim =  np.zeros((len(G), len(C), len(E), len(F)))
UpperLim =  np.zeros((len(G), len(C), len(E), len(F)))
arr_dS   =  np.zeros((len(G), len(C), len(E), len(F)))

LowerLim[:,:,:] = S0[:, :, :, :]
UpperLim[:,:,:] = S1[:, :, :, :]

UpperLim[UpperLim==SH]=np.nan
UpperLim[LowerLim<SH] =np.nan

LowerLim[LowerLim<SH]=np.nan
LowerLim[LowerLim==SH]=1

arr_dS[:,:,:,:]= SH*LowerLim[:,:,:,:]-UpperLim[:,:,:, :]

fig, ax = plt.subplots(ncols = 2, sharey=True)
#0
tmp_upper = np.ones((len(G),4))
tmp_lower = np.ones((len(G),3))
counter_plot=0
for ci in range(0,len(C)):
    for fi in range(0,len(F)):
        counter_plot+=1
        #1
        tmp_upper = np.ones((len(G),5))
        tmp_lower = np.ones((len(G),3))
        for gi in range(0,len(G)):
            data = np.squeeze(arr_dS[gi,ci,:,fi])
            data[np.isnan(data)]=-999
            tmp  = np.where(data!=-999)
            # 2 &3
            if tmp[0].size> 0: # then there are  entries to use
                ind= int(tmp[0][ 0])
                flx = (E[ind]*A0/(100*yr2sc))/Q[gi,ci,ind,fi]#10**(-6)*arrFA1[0,gi,ci,ind,fi]#/ (E[ind]*(100*yr2sc))
                clr = data[ind]
                FLX = (F[fi]*A0*E[ind]/(100*yr2sc))/FLUX[gi,ci, ind,fi]
                tmp_upper[gi,:] = [gi, ind, flx, clr, FLX]
                
                
                ind= int(tmp[0][ -1])
                flx = (E[ind]*A0/(100*yr2sc))/Q[gi,ci,ind,fi]#/ (E[ind]*(100*yr2sc))
                tmp_lower[gi,:] = [gi, ind, flx]
            else: 
                tmp_upper[gi,:] = [gi, -999, -999, -999, -999]
                tmp_lower[gi,:] = [gi, -999, -999]
        #4 slice and slay, but carefully
        tmp = np.where((tmp_upper[:,1]!=-999) & (tmp_upper[:,1]!=0))
        if tmp[0].size> 1:
            length_upper = np.where((tmp_upper[:,1]!=-999) & (tmp_upper[:,1]!=0))[0]
            length_lower = np.where((tmp_upper[:,1]!=-999) & (tmp_upper[:,1]!=0))[0]
            
            upper = np.array(tmp_upper[length_upper[0]:length_upper[-1]+1,:])
            lower = np.array(tmp_lower[length_lower[0]:length_lower[-1]+1,:])
            
            ind_u = np.array(tmp_upper[length_upper[0]:length_upper[-1]+1,1], dtype=int)
            ind_l = np.array(tmp_lower[length_lower[0]:length_lower[-1]+1,1], dtype=int)
           #5 fit this shit
            x_data = E[ind_u[:]]/100
            y_data = upper[:,2]
            #pre parameters
            a_tmp= (y_data[0]-y_data[-1])//((1/x_data[0]) - (1/x_data[-1]))
            c_tmp= y_data[-1] - a_tmp
            b_tmp= 0
            
            ax[0].scatter(x_data, y_data, c='k')

ax[0].set_xlim([0.6, 1.0])
ax[0].set_xticks([0.6, 0.7, 0.8, 0.9, 1.0])
ax[0].set_xticklabels([60, 70, 80, 90, 100])
ax[0].set_ylim([4,9])
ax[0].set_xlabel("net-evaporation [cm/yr]")
ax[0].set_ylabel("fwb / Q")

ax[0].set_yscale('linear' )
# ax[0].yaxis.set_major_locator(MultipleLocator(0.001))
# ax[0].yaxis.set_major_formatter('{x:.0f}')

ax[1].spines['right'].set_visible(False)
ax[1].spines['top' ].set_visible(False)
ax[1].spines['left'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)

ax[1].axes.get_xaxis().set_ticks([])

ax[0].spines['right'].set_visible(False)
ax[0].spines['top' ].set_visible(False)
ax[0].spines['left'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
fig.suptitle("Scenario A2 - raw")
plt.subplots_adjust(left=0.2,
                    bottom=0.1,
                    right=0.99,
                    top=0.75,
                    wspace=0.01,
                    hspace=0.4)
figurename= ("TMP_Scen_A2_ALL_new")
plt.savefig(figurename +"-raw.png", dpi=300)
plt.show()
plt.figure()

#%% Analysis
def model_f(x,a,b,c):
  return (a/(x-b)) + c
# Find Conditions when coeval precipitation would occur
# S0=Sh & S1<SH
LowerLim =  np.zeros((len(G), len(C), len(E), len(F)))
UpperLim =  np.zeros((len(G), len(C), len(E), len(F)))
arr_dS   =  np.zeros((len(G), len(C), len(E), len(F)))

LowerLim[:,:,:] = S0[:, :, :, :]
UpperLim[:,:,:] = S1[:, :, :, :]

UpperLim[UpperLim==SH]=np.nan
UpperLim[LowerLim<SH] =np.nan

LowerLim[LowerLim<SH]=np.nan
LowerLim[LowerLim==SH]=1

arr_dS[:,:,:,:]= SH*LowerLim[:,:,:,:]-UpperLim[:,:,:, :]

fig, ax = plt.subplots(ncols = 2, sharey=True)
#0
tmp_upper = np.ones((len(G),4))
tmp_lower = np.ones((len(G),3))
counter_plot=0
for ci in range(0,len(C)):
    for fi in range(0,len(F)):
        counter_plot+=1
        #1
        tmp_upper = np.ones((len(G),5))
        tmp_lower = np.ones((len(G),3))
        for gi in range(0,len(G)):
            data = np.squeeze(arr_dS[gi,ci,:,fi])
            data[np.isnan(data)]=-999
            tmp  = np.where(data!=-999)
            # 2 &3
            if tmp[0].size> 0: # then there are  entries to use
                ind= int(tmp[0][ 0])
                flx = (E[ind]*A0/(100*yr2sc))/Q[gi,ci,ind,fi]#10**(-6)*arrFA1[0,gi,ci,ind,fi]#/ (E[ind]*(100*yr2sc))
                clr = data[ind]
                FLX = (F[fi]*A0*E[ind]/(100*yr2sc))/FLUX[gi,ci, ind,fi]
                tmp_upper[gi,:] = [gi, ind, flx, clr, FLX]
                
                
                ind= int(tmp[0][ -1])
                flx = (E[ind]*A0/(100*yr2sc))/Q[gi,ci,ind,fi]#/ (E[ind]*(100*yr2sc))
                tmp_lower[gi,:] = [gi, ind, flx]
            else: 
                tmp_upper[gi,:] = [gi, -999, -999, -999, -999]
                tmp_lower[gi,:] = [gi, -999, -999]
        #4 slice and slay, but carefully
        tmp = np.where((tmp_upper[:,1]!=-999) & (tmp_upper[:,1]!=0))
        if tmp[0].size> 1:
            length_upper = np.where((tmp_upper[:,1]!=-999) & (tmp_upper[:,1]!=0))[0]
            length_lower = np.where((tmp_upper[:,1]!=-999) & (tmp_upper[:,1]!=0))[0]
            
            upper = np.array(tmp_upper[length_upper[0]:length_upper[-1]+1,:])
            lower = np.array(tmp_lower[length_lower[0]:length_lower[-1]+1,:])
            
            ind_u = np.array(tmp_upper[length_upper[0]:length_upper[-1]+1,1], dtype=int)
            ind_l = np.array(tmp_lower[length_lower[0]:length_lower[-1]+1,1], dtype=int)
           #5 fit this shit
            x_data = E[ind_u[:]]/100
            y_data = upper[:,2]
            #pre parameters
            a_tmp= (y_data[0]-y_data[-1])//((1/x_data[0]) - (1/x_data[-1]))
            c_tmp= y_data[-1] - a_tmp
            b_tmp= 0
            popt, pcov = curve_fit(model_f, x_data, y_data, p0=[a_tmp, b_tmp, c_tmp])
            a_opt, b_opt, c_opt = popt
            x_model = np.linspace(min(x_data), max(x_data), 500)
            y_model = model_f(x_model, a_opt, b_opt, c_opt) 
            
            f_Color     = interp1d( E[ind_u[:]]/100, upper[:,3], kind='slinear')
            
            pcm = ax[0].scatter(x_model, y_model, c= f_Color(x_model),
                                  norm=LogNorm(vmin=1, vmax=50))

            ax[1].annotate("c=10^" +str(int(np.log10(C[ci])))+
                           #"c= {:.2f}".format(C[ci]) +
                           ", {:.0f}".format(F[fi]*100) + "% of area \n-> (-fwb/F)*=10^["
                           + str(int(np.log10(min(upper[:,4]))))
                           + ", "+ str(int(np.log10(max(upper[:,4])))) + " ] ",
                           xy = (0,y_model[-200]), xycoords= 'data', 
                           xytext = (0.1,counter_plot/6), textcoords= 'axes fraction',
                           verticalalignment= 'top', horizontalalignment='left',
                           arrowprops=dict(facecolor='black', shrink=0.01, linewidth =0.00001))
            
fig.colorbar(pcm, ax= ax[:2], location= "top", 
           label = "salinity difference [kg/m³]")
ax[0].set_xlim([0.6, 1.0])
ax[0].set_xticks([0.6, 0.7, 0.8, 0.9, 1.0])
ax[0].set_xticklabels([60, 70, 80, 90, 100])
ax[0].set_ylim([4,9])
ax[0].set_xlabel("net-evaporation [cm/yr]")
ax[0].set_ylabel("fwb / Q")

ax[0].set_yscale('linear' )
# ax[0].yaxis.set_major_locator(MultipleLocator(0.001))
# ax[0].yaxis.set_major_formatter('{x:.0f}')

ax[1].spines['right'].set_visible(False)
ax[1].spines['top' ].set_visible(False)
ax[1].spines['left'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)

ax[1].axes.get_xaxis().set_ticks([])

ax[0].spines['right'].set_visible(False)
ax[0].spines['top' ].set_visible(False)
ax[0].spines['left'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
fig.suptitle("Scenario A2 - smoothened")
plt.subplots_adjust(left=0.2,
                    bottom=0.1,
                    right=0.99,
                    top=0.75,
                    wspace=0.01,
                    hspace=0.4)
figurename= ("TMP_Scen_A2_ALL_new")
plt.savefig(figurename +"-smooth.png", dpi=300)
plt.show()
plt.figure()
