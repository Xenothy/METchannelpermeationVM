#=============================================================================
#========================= INPUT THINGS HERE =================================
#=============================================================================

# This script takes four input files: the fractional block data and the time constant data
# for each of the two datasets to be compared. 
# Both should be CSV files, and should be located in the same folder as this script.
# Examples for how the CSVs should be formatted are in the folder. 
# They must be formatted exactly this way.


# The output is a 4-panel figure with every plot that describes the model, including:
# 1. Hill coefficients and Half-block concentrations, 2. Time constants, 
# 3. Concentration-dependent entry rates, 4. Energy profiles.


# N.B.: Panel 3. might look a little weird depending on the number of points 
# specified in lines 398-399 to create the curve. Increasing the number of points will
# make the curve smoother, but it will also take longer to run. 


# For running entry rate model
concentration_For_Model = 1
voltage_For_Model = -55 # in mV
number_Of_Channels = 1
p0 = 1 #resting open probability
T = 25  #temperature

alpha = 0.05 # for confidence intervals

# set the path to the folder on computer that contains the data
path = r"C:\Users\virgi\Documents\PhD\Experiments\Analysis\Data\Data for github" 

fixHillCo = 'y' #'y' for yes, 'n' to use the average calculated hillco across 8 
                # lowest voltages, 'fit' to use from fractional block fit
fixHillCoVal = 1


# Select which drug you want to compare apical and basal 

drug = 'gentamicin'
#drug = 'kanamycin'
#drug = 'amikacin'

figsize = 7,5
fontsize = 16
markersize = 8
linewidth = 2
capsize = 3
capthick = 1.5
elinewidth = 1.5
labels = 'Apical', 'Basal'

#============================================================================
# Import the packages

import numpy as np
import matplotlib.pylab as plt
import matplotlib.ticker as ticker
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import t

#=============================================================================
# Set up the plot

fig, axarr = plt.subplots(2,2, figsize = (5,5))
fig.tight_layout()

markers = ['o', 'v', 's','^', 'D', '<', 'P', '>', 'X']
colors = ['#FF0066', '#CC0099', '#770088', '#0000CC', '#00a0b4', '#00a033', '#d5b700', '#FF6600', '#CC0000']
markersize = 7


if drug == "gentamicin":
    kdylim = -10, 150
    nhylim = 2
    initial_halfBlock_guess = 10
    initial_nH_guess = 1
    initial_bindE_guess = -10
    color = ['#3399FF', 'b']
             
    tauYAxis = 0,7000
    tauXAxis = 0,3.2e-5
    tau_initial = 5E10

    vdfile1 = path + "\Gent_apical_frac.csv"
    vdfile2 = path + "\Gent_basal_frac.csv"
    tcfile1 = path + "\Gent_apical_tau.csv"
    tcfile2 = path + "\Gent_basal_tau.csv"
    
    
if drug == "kanamycin":
    kdylim = -500, 4100
    nhylim = 2
    initial_halfBlock_guess = 70
    initial_nH_guess = 1
    initial_bindE_guess = -5
    color = ['#22BB22', '#006400']
             
    tauYAxis = 0,60000
    tauXAxis = 0,11e-5
    tau_initial = 2.26e08
    
    vdfile1 = path + "\Kan_apical_frac.csv"
    vdfile2 = path + "\Kan_basal_frac.csv"
    tcfile1 = path + "\Kan_apical_tau.csv"
    tcfile2 = path + "\Kan_basal_tau.csv"
    
if drug == "amikacin":
    kdylim = -100, 2000
    nhylim = 2
    initial_halfBlock_guess = 20
    initial_nH_guess = 1
    initial_bindE_guess = -5
    color = ['#FF00FF', '#800080']
             
    tauYAxis = -1000,10000
    tauXAxis = 0,11e-5
    tau_initial = 1e8
             
             
    vdfile1 = path + "\Amik_apical_frac.csv"
    vdfile2 = path + "\Amik_basal_frac.csv"
    tcfile1 = path + "\Amik_apical_tau.csv"
    tcfile2 = path + "\Amik_basal_tau.csv"

#=============================================================================
#======================= Voltage dependence plot =============================
#=============================================================================


def doseResp(x, halfB, nH):
    term1 = 0
    term2 = (1 - 0)
    term3 = (1 + (x / halfB)** nH)
    return term1 + (term2 / term3)


axarr[0,0].xaxis.set_major_locator(ticker.MultipleLocator(20))
axarr[0,0].xaxis.set_minor_locator(ticker.MultipleLocator(10))
ax2 = axarr[0,0].twinx()


def DoVDPlot(file, color, num):
    frac = pd.read_csv(file, sep = ",", header = [0])
    frac.set_index(['Vm (mV)'], inplace =True)
    voltages = list(np.unique(frac.index.values))
    
    # Empty lists for hillCos from each loop.
    hillCo_list = []
    hillCo_CI_list = []
    halfBlock_list = []
    halfBlock_CI_list = []
    
    #===========================================================================
    # LOOP FOR DOSE RESPONSE CURVE FITTING
    # m is the voltage
    
    # Dose Response curve fitting function
    
    # Loop through each voltage
    for m in range (0, len(voltages)-5):
        xDoseRep = frac.loc[(voltages[m]), 'Concentration (uM)']
        yDoseRep = frac.loc[(voltages[m]), 'Current (pA)']
        
        initialDR = (initial_halfBlock_guess, initial_nH_guess)
        boundsDR = ((0, 0),
                    (1000000, 20))
     
        # Curve fitting each Dose Response
        pOpt, pCov = curve_fit(doseResp, xDoseRep, yDoseRep, p0 = initialDR, bounds = boundsDR, maxfev=10000)
        # Extracting hill coefficient, for later.
        hillCo_list.append(pOpt[1])
        halfBlock_list.append(pOpt[0])
        
        # Calculate confidence intervals
        hillCo_CI_list.append(np.sqrt(np.diag(pCov)[1]) * (t.ppf(1-alpha/2, df= len(xDoseRep)-len(pOpt))))
        halfBlock_CI_list.append(np.sqrt(np.diag(pCov)[0]) * (t.ppf(1-alpha/2, df= len(xDoseRep)-len(pOpt))))
        
    #=============================================================================
    # PLOT HALFBLOCK LIST
        
    axarr[0,0].errorbar(voltages[0:(m+1)], halfBlock_list, yerr = halfBlock_CI_list, linewidth = linewidth, capsize = capsize, capthick = capthick, ms=markersize, marker='o', color = color, label = labels[num], elinewidth = elinewidth)

    #=============================================================================
    # GET AND PLOT HILL COEFFICIENT
    
    ax2.errorbar(voltages[0:(m+1)], hillCo_list, yerr= hillCo_CI_list, linewidth = linewidth, capsize = capsize, capthick = capthick, marker='o', ms=markersize, mfc = 'w', color = color, elinewidth = elinewidth)

    hillCo = sum(hillCo_list) / len(hillCo_list)

    return voltages, halfBlock_list, halfBlock_CI_list, hillCo


apical_halfBlock_list = DoVDPlot(vdfile1, color[0], 0)[1]
apicalhb = pd.DataFrame(list(apical_halfBlock_list))
apicalhb = apicalhb.transpose()
apicalmaxhb = apicalhb[1].min()
apicalhillCo = DoVDPlot(vdfile1, color[0], 0)[3]

basal_halfBlock_list = DoVDPlot(vdfile2, color[1], 1)[1]
basalhb = pd.DataFrame(list(basal_halfBlock_list))
basalhb = basalhb.transpose()
basalmaxhb = basalhb[1].min()
basalhillCo = DoVDPlot(vdfile2, color[1], 0)[3]


axarr[0,0].set_ylabel(r"k$_{D}$ ($\mu$M) " + u" (\u26AB)")
axarr[0,0].set_xlabel("Voltage (mV)")
axarr[0,0].set_title("Voltage Dependence of Block")
ax2.set_ylabel(r"n$_{H}$" + u" (\u26AA)")

print ('Apical max half block =', apicalmaxhb)
print('\n')

print('\n')
print ('Basal max half block =', basalmaxhb)
print('\n')


#=============================================================================
#=========================== Fractional block fit ============================
#=============================================================================

def fracBlockFit(tupledat, hillCo, bindE, siteDist, zet, DE):
    (x,y) = tupledat
    term1 = (y*float(1.0e-9)*1000) ** hillCo
    term2 = np.exp(bindE + (siteDist * x * (1.0e-3 / 0.026) * zet))
    term3 = np.exp( -DE - (x * (1.0e-3 / 0.026) * zet))
    final = 1.0 / (1.0 + (term1 / (term2 * (1 + term3))))
    return final.ravel()

    
#prepare for the fit
if fixHillCo == 'y':
    initial = (fixHillCoVal, initial_bindE_guess, 1, 1.7, 5)
    parBounds = ((fixHillCoVal*0.99, -50, 0, 0, 0),
                (fixHillCoVal*1.01, 0, 5, 5, 30))
    
if fixHillCo == 'n':
    initial = (initial_nH_guess, initial_bindE_guess, 1, 1.7, 5)
    parBounds = ((0, -50, 0, 0, 0),
                (5, 0, 5, 5, 30))
    
if fixHillCo == 'fit':
    initial = (initial_nH_guess, initial_bindE_guess, 1, 1.7, 5)
    parBounds = ((0, -50, 0, 0, 0),
                (5, 0, 5, 5, 30))

def fracBlockDoThis(file, num):
    
    # re-import the csv but with concentration as the index.
    frac = pd.read_csv(file, sep = ",", header = [0])
    frac.set_index(['Concentration (uM)'], inplace =True)
    
    fitting_concs = frac.index.values
    fitting_volts = frac['Vm (mV)'].values
    fitting_currents = frac['Current (pA)'].values

    
    concentrations = list(np.unique(frac.index.values).astype("float"))
    voltages = list(np.unique(frac['Vm (mV)']).astype("float"))
    X, Y = np.meshgrid(fitting_volts, fitting_concs)
    numbSmoothDataPoints = 1000
    smoothX, smoothY = np.meshgrid(np.linspace(np.min(voltages), np.max(voltages), numbSmoothDataPoints),
                                   np.linspace(np.min(concentrations), np.max(concentrations), numbSmoothDataPoints))
    
    #do the fit!
    frac_pOpt, frac_pCov = curve_fit(fracBlockFit, (fitting_volts, fitting_concs), fitting_currents, p0 = initial, bounds = parBounds, maxfev=1000)
    
    return frac_pOpt



file1_frac = fracBlockDoThis(vdfile1, 0)
file2_frac = fracBlockDoThis(vdfile2, 1)

print("The Hill Coefficient was...\n")

#prepare for the fit
if fixHillCo == 'y':
    hillCo1 = fixHillCoVal
    hillCo2 = fixHillCoVal
    print("Fixed at ", fixHillCoVal)
if fixHillCo == 'fit':
    hillCo1 = file1_frac[0]
    hillCo2 = file2_frac[0]
    print("Fitted from fractional block curve: \n")
    print("    Apical: ,", hillCo1)
    print("    Basall: ,", hillCo2)
if fixHillCo == 'n':
    hillCo1 = apicalhillCo
    hillCo2 = basalhillCo
    print("Fitted from dose response curve list:  \n")
    print("    Apical: ,", hillCo1)
    print("    Basall: ,", hillCo2)
    
    
#=============================================================================
#========================= Time constants plot ===============================
#=============================================================================


def TimeConstantFit(x, m, b):
    return b + (m * x)

def TCPlotThis(file, color, hillCo):
    
    tcdata = pd.read_csv(file, sep = ",", header = [0])
    tauxData = ((tcdata['Concentration (uM)'])*(10**-6))**hillCo
    
    tauData = tcdata['1/Tau (Hz)']
    
    
    initialTC = (tau_initial, 0)
    boundsTC = ((-1E10, 0),
                (1.002E20, 100000))
    
    tau_x_grouped = tcdata.groupby('Concentration (uM)')['Concentration (uM)'].mean() #YES -fit
    tau_x_grouped = (tau_x_grouped*(10**-6))**hillCo #YES -fit
    
    avg_tau = tcdata.groupby('Concentration (uM)', as_index=False)['1/Tau (Hz)'].mean() #YES - fit

    
    tau_err = tcdata.groupby('Concentration (uM)', as_index=False)['Error (Hz)']
    tau_CI = tau_err.sem()*(t.ppf(1-alpha/2, df=(len(avg_tau)-1)))
    
    ### FOR FITTING TO THE RAW DATA:
    raw_tau_pOpt, raw_tau_pCov = curve_fit(TimeConstantFit, tauxData, tauData, p0 = initialTC, bounds = boundsTC, maxfev = 100000)
    

    tau_slope = raw_tau_pOpt[0]
    tau_slope_CI = np.sqrt((np.diag(raw_tau_pCov)[0])) * (t.ppf(1-alpha/2, df= (len(tauData)-len(raw_tau_pOpt))))
    
    
    xvals = np.linspace(np.min(tauxData), np.max(tauxData), 1000)
    fitYVals = TimeConstantFit(xvals, raw_tau_pOpt[0], raw_tau_pOpt[1])
    
    tc_x = ((avg_tau.iloc[:,0])*(10**-6))**hillCo
    tc_y = avg_tau.iloc[:,1]
    tc_err = tau_CI.iloc[:,1]
    
    axarr[0,1].errorbar(tc_x, tc_y, yerr = tc_err, fmt='o', marker='o', color=color, capsize=capsize, capthick=capthick, elinewidth = elinewidth, ms = markersize)
    axarr[0,1].plot(xvals, fitYVals, color = color, lw = linewidth)

    
    return tau_slope, tau_slope_CI


TCPlotThis(tcfile1, color[0], hillCo1)
TCPlotThis(tcfile2, color[1], hillCo2)

axarr[0,1].set_xlabel("Concentration (M) ^ nH")
axarr[0,1].set_ylabel('1 / Tau (Hz)')
axarr[0,1].set_title("Time Constants")
axarr[0,1].set_ylim(tauYAxis)
axarr[0,1].set_xlim(tauXAxis)

# Constants
kT = (273.15 + T ) * 1.38e-23 #(Boltzman's constant) , absolute temperature
k0 = 6.18401E+12  # kT/6.63e-34, thermal noise energy divided Planck's Constant, proportionality constant of the first order rate constant. 

# fitted values
Eb1 = file1_frac[1]  # In terms of kT binding site
Deltab1 = file1_frac[2] # binding site relative distance from first energy barrier
DeltaE1 = file1_frac[4] # difference in energy barriers
zet1 = file1_frac[3] # charge of the molecule

Eb2 = file2_frac[1]  # In terms of kT binding site
Deltab2 = file2_frac[2] # binding site relative distance from first energy barrier
DeltaE2 = file2_frac[4] # difference in energy barriers
zet2 = file2_frac[3] # charge of the molecule

#from slope of time constant fit (in per mole per second ^ nhill)
k11 = TCPlotThis(tcfile1, color[0], hillCo1)[0]
k12 = TCPlotThis(tcfile2, color[1], hillCo2)[0]

print('The apical time constant is:', k11)
print('The basal time constant is:', k12)

print('\n\n')

Nch = number_Of_Channels #number of channels
D = concentration_For_Model  # in mM

def MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, uM, mV):  
    V = mV * (10**-3)
    D = uM * (10**3)
    Vs = kT / zet / 1.6e-19
    K1 = np.exp(Eb + Deltab * V / Vs) * (1 + np.exp(-DeltaE - V / Vs))
    E1 = -np.log(k1 / k0)  # In terms of kT
    E2 = E1 + DeltaE
    k2 = k0 * np.exp(-(E2 - Eb) - (1 - Deltab) * V / Vs)
    Nentry = Nch * p0 * k2 / (1 + K1 / (D * 1e-9) ** n) # per cell
    return Nentry, E1, E2

#Nentry_concentrations = np.linspace(0.01, 100000, 1000000) # for best results
Nentry_concentrations = np.linspace(0.01, 100000, 10000) # to speed up the modelling

Nentry_conc_list1 = []
Nentry_conc_list2 = []

for p in range (0, len(Nentry_concentrations)):
    Nentry_conc1 = MET_Model(Eb1, Deltab1, DeltaE1, zet1, hillCo1, k11, p0, T, kT, k0, Nch, Nentry_concentrations[p], voltage_For_Model)[0]
    Nentry_conc_list1.append(Nentry_conc1)
    
    Nentry_conc2 = MET_Model(Eb2, Deltab2, DeltaE2, zet2, hillCo2, k12, p0, T, kT, k0, Nch, Nentry_concentrations[p], voltage_For_Model)[0]
    Nentry_conc_list2.append(Nentry_conc2)

axarr[1,0].plot(Nentry_concentrations, Nentry_conc_list1, color = color[0])
axarr[1,0].plot(Nentry_concentrations, Nentry_conc_list2, color = color[1])
axarr[1,0].set_xlabel('Concentration ($\mu$M)')
axarr[1,0].set_ylabel('Entry Rate (molecules/s)')
axarr[1,0].semilogx()
axarr[1,0].set_title('Entry Rates')

#=============================================================================
#========================= Energy profiles plot ==============================
#=============================================================================

def EnergyProfile(x, E1, d1, w1, Eb, db, wb, E2, d2, w2):
    term1 = 4*E1*1/(1+np.exp(-(x-d1)/w1))*(1-1/(1+np.exp(-(x-d1)/w1)))
    term2 = 4*Eb*1/(1+np.exp(-(x-db)/wb))*(1-1/(1+np.exp(-(x-db)/wb)))
    term3 = 4*E2*1/(1+np.exp(-(x-d2)/w2))*(1-1/(1+np.exp(-(x-d2)/w2)))
    term1 + term2 + term3
    return term1 + term2 + term3

E11 = MET_Model(Eb1, Deltab1, DeltaE1, zet1, hillCo1, k11, p0, T, kT, k0, Nch, concentration_For_Model, voltage_For_Model)[1]
E21 = MET_Model(Eb1, Deltab1, DeltaE1, zet1, hillCo1, k11, p0, T, kT, k0, Nch, concentration_For_Model, voltage_For_Model)[2]

E12 = MET_Model(Eb2, Deltab2, DeltaE2, zet2, hillCo2, k12, p0, T, kT, k0, Nch, concentration_For_Model, voltage_For_Model)[1]
E22 = MET_Model(Eb2, Deltab2, DeltaE2, zet2, hillCo2, k12, p0, T, kT, k0, Nch, concentration_For_Model, voltage_For_Model)[2]

xRange = np.linspace(-0.1, 1.1, 1000)

apical_ep = EnergyProfile(xRange, E11, 0, 0.01, Eb1, Deltab1, 0.025, E21, 1, 0.01)
basal_ep = EnergyProfile(xRange, E12, 0, 0.01, Eb2, Deltab2, 0.025, E22, 1, 0.01)

axarr[1,1].plot(xRange, apical_ep, color = color[0], lw = linewidth)
axarr[1,1].plot(xRange, basal_ep, color = color[1], lw = linewidth)
axarr[1,1].set_ylabel('Free Energy Difference (kT)')
axarr[1,1].set_xlabel('Relative Electrical Distance')
axarr[1,1].set_title('Energy Profiles')
axarr[1,1].legend(labels, frameon = False)

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()