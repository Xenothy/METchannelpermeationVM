#=============================================================================
#========================= INPUT THINGS HERE =================================
#=============================================================================

# This script takes two input files: the fractional block data and the time constant data. 
# Both should be CSV files, and should be located in the same folder as this script.
# Examples for how the CSVs should be formatted are in the folder. 
# They must be formatted exactly this way.


# The output is an 8-panel figure with every plot that describes the model, including:
# 1. Dose response curves, 2. Hill coefficients, 3. Half-block concentrations, 
# 4. Time constants, 5. Fractional block curves, 6. Concentration-dependent entry rate, 
# 7. Energy profile, 8. 3D  entry rate

# N.B.: Panel 6. might look a little weird depending on the number of points 
# specified in lines 550-551 to create the curve. Increasing the number of points will
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

# Setting the hill coefficient
fixHillCo = 'y' #'y' for yes, 'n' to use the average calculated hillco across 8 
                # lowest voltages, 'fit' to use from fractional block fit
fixHillCoVal = 1


# Select which drug
drug = 'gentamicin'
#drug = 'kanamycin'
#drug = 'amikacin'


# Select which coil
coil = 'apical'
#coil = 'basal'
#
#============================================================================
# Import the packages

import numpy as np
import matplotlib.pylab as plt
import matplotlib.ticker as ticker
from matplotlib import container
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import t
from matplotlib.ticker import LogFormatter
from mpl_toolkits.mplot3d import Axes3D # yes it is used ignore the error

#=============================================================================
# Set up the plot

#fig, axarr = plt.subplots(2,4, gridspec_kw={'width_ratios': [3,2,3,3]})
fig, axarr = plt.subplots(2,4, figsize = (15,5))
fig.tight_layout()

markers = ['o', 'v', 's','^', 'D', '<', 'P', '>', 'X']
colors = ['#FF0066', '#CC0099', '#770088', '#0000CC', '#00a0b4', '#00a033',
          '#d5b700', '#FF6600', '#CC0000']
markersize = 7

#============================================================================
#====================== Setting up the conditions ===========================
#============================================================================

if drug == "gentamicin":
    initial_bindE_guess = -10
    initial_halfBlock_guess = 10
    initial_nH_guess = 1
    tauYAxis = 0,7000
    tauXAxis = 0,3.2e-5
    tau_initial = 1e8
    hillCo_init = 1
    
    if coil == "apical":
            file = path + "\Gent_apical_frac.csv"
            file2 = path + "\Gent_apical_tau.csv"
            color = '#3399FF'
            
    elif coil == "basal":
            file = path + "\Gent_basal_frac.csv"
            file2 = path + "\Gent_basal_tau.csv"
            color = 'b'
            
if drug == "kanamycin":
    initial_bindE_guess = -5
    initial_halfBlock_guess = 70
    initial_nH_guess = 1
    tauYAxis = 0,60000
    tauXAxis = 0,11e-5
    tau_initial = 2.26e08
    hillCo_init = 1
    
    if coil == "apical":
            file = path + "\Kan_apical_frac.csv"
            file2 = path + "\Kan_apical_tau.csv"
            color = '#22BB22'
            
    elif coil == "basal":
            file = path + "\Kan_basal_frac.csv"
            file2 = path + "\Kan_basal_tau.csv"
            color = '#006400'
            
if drug == "amikacin":
    initial_bindE_guess = -5
    initial_halfBlock_guess = 20
    initial_nH_guess = 1
    tauYAxis = -1000,10000
    tauXAxis = 0,11e-5
    tau_initial = 1e8
    hillCo_init = 1
             
    if coil == "apical":
            file = path + "\Amik_apical_frac.csv"
            file2 = path + "\Amik_apical_tau.csv"
            color = '#FF00FF'
            
    elif coil == "basal":
            file = path + "\Amik_basal_frac.csv"
            file2 = path + "\Amik_basal_tau.csv"
            color = '#800080'


print(file)
print(file2)
print("\n")


#=============================================================================
#=================== Dose Response Curve Fitting =============================
#=============================================================================


# import frac CSV, with header as first row
frac = pd.read_csv(file, sep = ",", header = [0])
frac.set_index(['Vm (mV)'], inplace =True)

## get a list of the individual concentrations and voltages
concentrations = list(np.unique(frac['Concentration (uM)']))
voltages = list(np.unique(frac.index.values))

# Empty lists to get hillCos and halfblocks from each loop.
hillCo_list = []
hillCo_CI_list = []
halfBlock_list = []
halfBlock_CI_list = []

# Dose Response curve fitting function

# Here, the upper and lower bound are fixed to 1 and 0, respectively. If you 
# want to change that, change 1 and and 0 to parameters and add bounds and 
# initialising paramters below.

def doseResp(x, halfB, nH):
    term1 = 0
    term2 = (1 - 0)
    term3 = (1 + (x / halfB)** nH)
    return term1 + (term2 / term3)

# Loop through each voltage
for m in range (0, len(voltages)-5):
    
    xDoseRep = frac.loc[(voltages[m]), 'Concentration (uM)']
    yDoseRep = frac.loc[(voltages[m]), 'Current (pA)']

    initialDR = (initial_halfBlock_guess, initial_nH_guess)
    boundsDR = ((0, 0),
                (30000, 10))
 
    # Curve fitting each Dose Response
    pOpt, pCov = curve_fit(doseResp, xDoseRep, yDoseRep, p0 = initialDR, 
                           bounds = boundsDR, maxfev=10000)
    
    hillCo_list.append(pOpt[1])
    halfBlock_list.append(pOpt[0])
    
    # Calculate confidence intervals
    hillCo_CI_list.append(np.sqrt(np.diag(pCov)[1]) * (t.ppf(1-alpha/2, 
                          df= len(yDoseRep)-len(pOpt))))
    halfBlock_CI_list.append(np.sqrt(np.diag(pCov)[0]) * (t.ppf(1-alpha/2, 
                             df= len(yDoseRep)-len(pOpt))))
    
    # Plotting with a nice smoothed function
    xRange = np.linspace(np.min(concentrations), np.max(concentrations), 1000) 
    axarr[0,0].plot(xRange, doseResp(xRange, *pOpt), color = colors[m]) 

print("The dose response hill coefficient is {:3.4}".format(sum(hillCo_list) / len(hillCo_list)))

    
# For plotting mean data
grouped_dr = frac.groupby(['Vm (mV)', 'Concentration (uM)'], as_index=True)['Current (pA)']
grouped_dr_mean = grouped_dr.mean()
grouped_dr_CI = grouped_dr.sem()*(t.ppf(1-alpha/2, df=(len(grouped_dr)-1)))
grouped_dr_n = grouped_dr.count()

dr_indices = np.unique(grouped_dr_mean.index.get_level_values('Vm (mV)'))

for j, i in enumerate(dr_indices[0:9]):
    if j != len(concentrations)-1: 
        axarr[0,0].errorbar(grouped_dr_mean[i].index, grouped_dr_mean[i].values, 
             yerr=grouped_dr_CI[i], color = colors[j], marker= markers[j], 
             label = "{} mV".format(i), capsize = 2, capthick = 2, ls ='none', 
             ms = markersize, antialiased = True)

# for plotting the fit - this removes error bars from the legend
handles, labels = axarr[0,0].get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]

axarr[0,0].set_xlabel("Concentration ($\mu$M)")
axarr[0,0].set_ylabel(r"I $_{drug}$ / I $_{control}$")
axarr[0,0].set_title('Dose Responses')
axarr[0,0].legend(handles, labels, loc = 'best', frameon=False, labelspacing = 0)
axarr[0,0].set_xlim((np.min(concentrations)-(0.5*np.min(concentrations))), 
     (np.max(concentrations)+ (0.5*np.max(concentrations))))
axarr[0,0].semilogx()


#=============================================================================
#===================== Fractional block curve fitting ========================
#=============================================================================

def fracBlockFit(tupledat, hillCo, bindE, siteDist, zet, DE):
    (x,y) = tupledat
    term1 = (y*float(1.0e-9)*1000) ** hillCo
    term2 = np.exp(bindE + (siteDist * x * (1.0e-3 / 0.026) * zet))
    term3 = np.exp( -DE - (x * (1.0e-3 / 0.026) * zet))
    final = 1.0 / (1.0 + (term1 / (term2 * (1 + term3))))
    return final.ravel()

# re-import the csv but with concentration as the index.
frac = pd.read_csv(file, sep = ",", header = [0])
frac.set_index(['Concentration (uM)'], inplace =True)

fitting_concs = frac.index.values
fitting_volts = frac['Vm (mV)'].values
fitting_currents = frac['Current (pA)'].values

#prepare for the fit
if fixHillCo == 'y':
    initial = (fixHillCoVal, initial_bindE_guess, 1, 1.7, 5)
    parBounds = ((fixHillCoVal*0.99, -50, 0, 0, 0),
                (fixHillCoVal*1.01, 0, 5, 5, 30))
    
if fixHillCo == 'n':
    initial = (hillCo_init, initial_bindE_guess, 1, 1.7, 5)
    parBounds = ((0, -50, 0, 0, 0),
                (5, 0, 5, 5, 30))
    
if fixHillCo == 'fit':
    initial = (hillCo_init, initial_bindE_guess, 1, 1.7, 5)
    parBounds = ((0, -50, 0, 0, 0),
                (5, 0, 5, 5, 30))

concentrations = list(np.unique(frac.index.values).astype("float"))
voltages = list(np.unique(frac['Vm (mV)']).astype("float"))
X, Y = np.meshgrid(fitting_volts, fitting_concs)
numbSmoothDataPoints = 1000
smoothX, smoothY = np.meshgrid(np.linspace(np.min(voltages), np.max(voltages), 
                                           numbSmoothDataPoints), np.linspace(np.min(concentrations), 
                                                                             np.max(concentrations), numbSmoothDataPoints))

#do the fit!
frac_pOpt, frac_pCov = curve_fit(fracBlockFit, (fitting_volts, fitting_concs), 
                                 fitting_currents, p0 = initial, bounds = parBounds, maxfev=1000)

frac_CI = (np.sqrt(np.diag(frac_pCov)))*(t.ppf(1-alpha/2, df=(len(fitting_currents)-1)))

print("\n\nFractional block fit results: \n")
print("The hill coefficiant fitted from fractional block curves is {:3.4},".format(frac_pOpt[0]), "with CI: {:3.4},".format(frac_CI[0]))
print("The binding energy is {:3.4},".format(frac_pOpt[1]), "with CI: {:3.4},".format(frac_CI[1]))
print("The site distance is {:3.4},".format(frac_pOpt[2]), "CI: {:3.4},".format(frac_CI[2]))
print("The charge of the molecule is {:3.4},".format(frac_pOpt[3]), "CI: {:3.4},".format(frac_CI[3]))
print("The difference in energy barriers is {:3.4},".format(frac_pOpt[4]), "CI: {:3.4},".format(frac_CI[4]))


#=============================================================================
#================== Fractional block plotting ================================
#=============================================================================

# to plot the raw data
#for i in range(0, len(concentrations)):
#    appendData = frac.loc[(concentrations[i]), 'Current (pA)'].astype("float").tolist()
#    xFracData = frac.loc[(concentrations[i]), 'Vm (mV)']
#    axarr[1,0].scatter(xFracData, appendData, marker='o')


#to plot the means of the data
frac_grouped = frac.groupby([ 'Concentration (uM)','Vm (mV)'], as_index=True)['Current (pA)']
frac_grouped_mean = frac_grouped.mean()
frac_grouped_CI = (frac_grouped.sem())*(t.ppf(1-alpha/2, df=(len(frac_grouped)-1)))
frac_grouped_n = frac_grouped.count()

frac_indices = np.unique(frac_grouped_mean.index.get_level_values('Concentration (uM)'))
colors.reverse()

for j, i in enumerate(frac_indices[0:len(concentrations)]):
    axarr[1,0].errorbar(frac_grouped_mean[i].index, frac_grouped_mean[i].values, yerr=frac_grouped_CI[i], 
            color = colors[j], label = r"{} $\mu$M (n={})".format(i, int(frac_grouped_n[(j*14):((j*14)+13)].mean())), 
            marker= markers[j], ms=markersize, capsize = 2, capthick = 2, ls = 'none', antialiased = True)
    

#to plot the smooth fitted line on top
def fracBlockPlot(x, conc, hillCo, bindE, siteDist, zet, DE):
    term1 = (conc*1.0e-9*1000) ** hillCo
    term2 = np.exp(bindE + (siteDist * x * (1.0e-3 / 0.026) * zet))
    term3 = np.exp( -DE - (x * (1.0e-3 / 0.026) * zet))
    return  1.0 / (1.0 + (term1 / (term2 * (1 + term3))))

voltRange = np.linspace(np.min(voltages), np.max(voltages), 1000)
plottingData_list = []
for p in range (0, len(concentrations)):
    plottingData_list = []
    for l in range (0, len(voltRange)):
        plottingData = fracBlockPlot(voltRange[l], concentrations[p], *frac_pOpt)
        plottingData_list.append(plottingData)
    axarr[1,0].plot(voltRange, plottingData_list, color = colors[p], lw = 0.75, antialiased = True)

# again removing error bars from the legend
handles, labels = axarr[1,0].get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]

# formatting the plot
axarr[1,0].xaxis.set_major_locator(ticker.MultipleLocator(40))
axarr[1,0].xaxis.set_minor_locator(ticker.MultipleLocator(20))
axarr[1,0].legend(handles, labels, loc = 'best', frameon=False,  fontsize = 10, labelspacing = 0)
axarr[1,0].set_xlabel('Voltage (mV)')
axarr[1,0].set_ylabel('Current (pA)')
axarr[1,0].set_ylim((0, 1.2))
axarr[1,0].set_title('Fractional Block')


#=============================================================================
#================= Get and plot hill coefficient =============================
#=============================================================================

print("\n")
print("The hill coefficient used in the model is:")

if fixHillCo == 'y':
    hillCo = fixHillCoVal
    print("Fixed at ", fixHillCoVal)
if fixHillCo == 'fit':
    hillCo = frac_pOpt[0]
    print("Fitted from fractional block curve: ", hillCo)
if fixHillCo == 'n':
    hillCo = sum(hillCo_list) / len(hillCo_list)
    print("Fitted from dose response curve list: ", hillCo)

axarr[0,1].errorbar(voltages[0:(m+1)], hillCo_list, yerr= hillCo_CI_list, 
     capsize=3, marker='o', ms=markersize, mfc = 'w', color = color)
axarr[0,1].set_ylim(0, hillCo + (hillCo*.5))
axarr[0,1].set_xlabel("Voltage (mV)")
axarr[0,1].set_ylabel(r"n$_{H}$")
axarr[0,1].set_title('Hill Coefficients')



#=============================================================================
#============= Plot half-blocking concentration list ========================
#=============================================================================


axarr[0,2].errorbar(voltages[0:(m+1)], halfBlock_list, yerr = halfBlock_CI_list, 
     capsize=3, ms=markersize, marker='o', color = color)
axarr[0,2].set_xlabel("Voltage (mV)")
axarr[0,2].set_ylabel(r"k$_{D}$ ($\mu$M)")
axarr[0,2].set_title('Half-Blocking Concentration')

#=============================================================================
#========================== Time Constants ===================================
#=============================================================================

def TimeConstantFit(x, m, b):
    return b + (m * x)

tcdata = pd.read_csv(file2, sep = ",", header = [0])
tauxData = ((tcdata['Concentration (uM)'])*(10**-6))**hillCo

initialTC = (tau_initial, 0)
boundsTC = ((-1e10, 0),
            (1.002E20, 100000))

#=============================================================================

xvals = np.linspace(np.min(tauxData), np.max(tauxData), 1000)**hillCo

tauX = tcdata['Concentration (uM)']*(10**-6)
tauData = tcdata['1/Tau (Hz)']
tauError = tcdata['Error (Hz)']

#axarr[0,3].errorbar(tauX, tauData, yerr = tauError, fmt='o', marker='o', capsize=5, capthick=3, color = color)

#============= FOR FITTING TO THE RAW DATA:

#pOpt_raw, tau_pCov_raw = curve_fit(TimeConstantFit, tauxData, tauData, p0 = initialTC, sigma= tauError,
#                           bounds = boundsTC, maxfev = 100000)


#============= FOR FITTING TO THE RAW, NO ERRORS:

pOpt_raw, pCov_raw = curve_fit(TimeConstantFit, tauxData, tauData, p0 = initialTC,
                           bounds = boundsTC, maxfev = 100000)

tauDataraw = tcdata['1/Tau (Hz)']
tau_raw_CI = (np.sqrt(np.diag(pCov_raw))) * t.ppf(1-alpha/2, df= (len(tauDataraw)-len(pOpt_raw)))

#============= FOR FITTING TO THE MEANS OF THE DATA:
tau_x_grouped = tcdata.groupby('Concentration (uM)')['Concentration (uM)'].mean()
tau_x_grouped = (tau_x_grouped*(10**-6))**hillCo
avg_tau = tcdata.groupby('Concentration (uM)', as_index=False)['1/Tau (Hz)'].mean()
tau_err = tcdata.groupby('Concentration (uM)', as_index=False)['Error (Hz)']

tcdata['ErrSq'] =(tcdata['Error (Hz)'])**2

tau_err_prop = np.sqrt(tcdata.groupby('Concentration (uM)', as_index=False)['ErrSq'].sum())
error_tau = tau_err_prop['ErrSq']

avg_tau_list = avg_tau['1/Tau (Hz)'].tolist()
avg_tau_error_list = error_tau.values

pOpt, tau_pCov = curve_fit(TimeConstantFit, tau_x_grouped, avg_tau_list, 
                       p0 = initialTC, bounds = boundsTC, 
                       sigma = avg_tau_error_list, maxfev = 100000)


#=============================================================================


tauSlope = pOpt_raw[0] # to use the slope from raw TC data
#tauSlope = pOpt_mean[0] # to use the slope time constant from mean TC data

tau_slope_CI = np.sqrt((np.diag(pCov_raw)[0])) * (t.ppf(1-alpha/2, df= (len(tauData)-len(pOpt_raw))))

tau_CI = tau_err.sem()*(t.ppf(1-alpha/2, df=(len(avg_tau)-1)))

tc_x = ((avg_tau.iloc[:,0])*(10**-6))**hillCo
tc_y = avg_tau.iloc[:,1]
tc_err = tau_CI.iloc[:,1]


fitYVals = TimeConstantFit(xvals, pOpt[0], pOpt[1])
 

axarr[0,3].errorbar(tc_x, tc_y, yerr = tc_err, fmt='o', marker='o', color= color, 
     ms=markersize, capsize = 2, capthick = 2, ls = 'none', antialiased = True)

axarr[0,3].plot(xvals, fitYVals, color = color)


axarr[0,3].set_xlabel('Concentration (M) ^ nH')
axarr[0,3].set_ylabel('1 / Tau (Hz)')
axarr[0,3].set_ylim(tauYAxis)
axarr[0,3].set_xlim(tauXAxis)
axarr[0,3].set_title('Time Constants')
axarr[0,3].ticklabel_format(style='sci', axis='x', scilimits=(0,0))

print("\n\n")
print('The time constant slope (k1) is {:3.4}'.format(tauSlope), 'with CI {:3.4}'.format(tau_slope_CI))

#=============================================================================
#============================= Modelling =====================================
#=============================================================================

# Constants
kT = (273.15 + T ) * 1.38e-23 #(Boltzman's constant) , absolute temperature
k0 = 6.18401E+12  # kT/6.63e-34, thermal noise energy divided Planck's Constant, proportionality constant of the first order rate constant. 

# fitted values
Eb = frac_pOpt[1]  # In terms of kT binding site
Eb_CI = frac_CI[1]
Deltab = frac_pOpt[2] # binding site relative distance from first energy barrier
Deltab_CI = frac_CI[2]
DeltaE = frac_pOpt[4] # difference in energy barriers
DeltaE_CI = frac_CI[4]
zet = frac_pOpt[3] # charge of the molecule
zet_CI = frac_CI[3]

#n = frac_pOpt[0] # hill coefficient
n = hillCo

k1 = tauSlope #from slope of time constant fit (in per mole per second ^ nhill)

k1_CI = tau_slope_CI
Nch = number_Of_Channels #number of channels
D = concentration_For_Model  # in nM, picking a concentration for now, could loop in future. 

# model!
def MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, uM, mV):  
    V = mV * (10**-3)
    D = uM * (10**3)
    Vs = kT / zet / 1.6e-19
    V0 = -Vs * (DeltaE + np.log(Deltab / (1 - Deltab)))
    K1 = np.exp(Eb + Deltab * V / Vs) * (1 + np.exp(-DeltaE - V / Vs))
    E1 = -np.log(k1 / k0)  # In terms of kT
    E2 = E1 + DeltaE
    k_1 = k0 * np.exp(-(E1 - Eb) + Deltab * V / Vs)
    k2 = k0 * np.exp(-(E2 - Eb) - (1 - Deltab) * V / Vs)
    k_2 = k0 * np.exp(-E2)
    k2 = k0 * np.exp(-((E1 + DeltaE) - Eb) - (1 - Deltab) * V / Vs)
    Nentry = Nch * p0 * k2 / (1 + K1 / (D * 1e-9) ** n)
    Half_Nentry = (Nch * p0 *k2)/2
    return Nentry, V0, E1, E2, k_1, k_2, k2, Half_Nentry

Nentry1 = MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, 1, voltage_For_Model)[0]
Nentry2 = MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, 2, voltage_For_Model)[0]
Nentry3 = MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, 3, voltage_For_Model)[0]
Nentry100 = MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, 100, voltage_For_Model)[0]
Nentry10000 = MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, 10000, voltage_For_Model)[0]


V0 = MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, D, voltage_For_Model)[1]
V0_CI = abs(V0 - MET_Model(Eb- Eb_CI, Deltab - Deltab_CI, DeltaE - DeltaE_CI, zet - zet_CI, n, k1 - k1_CI, p0, T, kT, k0, Nch, D, voltage_For_Model)[1])

E1 = MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, D, voltage_For_Model)[2]
E1_CI = abs(E1 - MET_Model(Eb- Eb_CI, Deltab - Deltab_CI, DeltaE - DeltaE_CI, zet - zet_CI, n, k1 - k1_CI, p0, T, kT, k0, Nch, D, voltage_For_Model)[2])

E2 = MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, D, voltage_For_Model)[3]
E2_CI = abs(E2 - MET_Model(Eb- Eb_CI, Deltab - Deltab_CI, DeltaE - DeltaE_CI, zet - zet_CI, n, k1 - k1_CI, p0, T, kT, k0, Nch, D, voltage_For_Model)[3])

k_1 = MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, D, voltage_For_Model)[4]
k_1_CI = abs(k_1 - MET_Model(Eb- Eb_CI, Deltab - Deltab_CI, DeltaE - DeltaE_CI, zet - zet_CI, n, k1 - k1_CI, p0, T, kT, k0, Nch, D, voltage_For_Model)[4])

k_2 = MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, D, voltage_For_Model)[5]
k_2_CI = abs(k_2 - MET_Model(Eb- Eb_CI, Deltab - Deltab_CI, DeltaE - DeltaE_CI, zet - zet_CI, n, k1 - k1_CI, p0, T, kT, k0, Nch, D, voltage_For_Model)[5])

k2 = MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, D, voltage_For_Model)[6]
k2_CI = abs(k2 - MET_Model(Eb- Eb_CI, Deltab - Deltab_CI, DeltaE - DeltaE_CI, zet - zet_CI, n, k1 - k1_CI, p0, T, kT, k0, Nch, D, voltage_For_Model)[6])


#=============================================================================
#=============== Plot concentration-dependent entry ==========================
#=============================================================================

# the number of points in the curve (3rd number in the two next lines) 
# will affect the maximum and half-entry rate calculations as well as
# how the entry rate curve looks. The more points, the more accurate, 
# but it will also make it run slower.
# Use fewer points to speed up the calculation, or more to make it look nicer.
 
#Nentry_concentrations = np.linspace(0.01, 100000, 1000000) # for best results
Nentry_concentrations = np.linspace(0.01, 100000, 10000) # to speed up the model

Nentry_conc_list = []

for p in range (0, len(Nentry_concentrations)):
    Nentry_conc = MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, Nentry_concentrations[p], voltage_For_Model)[0]
    Nentry_conc_list.append(Nentry_conc)
    
Max_entry = max(Nentry_conc_list)
Half_entry = round((max(Nentry_conc_list)/2), 2)


def find_nearest(Nentry_conc_list, number):
    smallestDif = 99999999
    smallestDifPos = 0
    for i in range(0, len(Nentry_conc_list)):
        dif = abs(Nentry_conc_list[i] - number) 
        if dif < smallestDif:
            smallestDif = dif
            smallestDifPos = Nentry_concentrations[i]
    return smallestDifPos


Half_entry_concentration = find_nearest(Nentry_conc_list, max(Nentry_conc_list)/2)
Max_entry_concentration = find_nearest(Nentry_conc_list, 0.99*max(Nentry_conc_list))

axarr[1,1].plot(Nentry_concentrations, Nentry_conc_list, color = color)
axarr[1,1].set_xlabel('Concentration ($\mu$M)')
axarr[1,1].set_ylabel('Entry Rate (molecules/s)')
axarr[1,1].semilogx()
axarr[1,1].set_title('Entry Rate')


#=============================================================================
#========================= Energy profile plot ===============================
#=============================================================================

def EnergyProfile(x, E1, d1, w1, Eb, db, wb, E2, d2, w2):
    term1 = 4*E1*1/(1+np.exp(-(x-d1)/w1))*(1-1/(1+np.exp(-(x-d1)/w1)))
    term2 = 4*Eb*1/(1+np.exp(-(x-db)/wb))*(1-1/(1+np.exp(-(x-db)/wb)))
    term3 = 4*E2*1/(1+np.exp(-(x-d2)/w2))*(1-1/(1+np.exp(-(x-d2)/w2)))
    term1 + term2 + term3
    return term1 + term2 + term3

epRange = np.linspace(-0.1, 1.1, 1000)
energy_profile = EnergyProfile(epRange, E1, 0, 0.01, Eb, Deltab, 0.01, E2, 1, 0.01)
axarr[1,2].plot(epRange, energy_profile, color = color)
axarr[1,2].set_ylabel('Free Energy difference (kT)')
axarr[1,2].set_xlabel('Relative electrical distance')


#=============================================================================
#========================= 3D entry rate plot ================================
#=============================================================================

concentrations = np.linspace(1, 1000, 1000)
voltages = np.linspace(-200, 0, 1000)
X,Y = np.meshgrid(voltages, concentrations)

Nentry = MET_Model(Eb, Deltab, DeltaE, zet, n, k1, p0, T, kT, k0, Nch, Y, X)[0]

ax = fig.add_subplot(2, 4, 8, projection='3d')

def makeLogAxesLabels(minTick, maxTick):
    listt = []
    for i in range(minTick, maxTick):
        temp = str((10**i)/10)
        listt.append(temp)
    return listt

maxYOrderOfTen = 5

myLogAxesList = makeLogAxesLabels(0, maxYOrderOfTen)

ax.plot_wireframe(X, np.log10(Y), Nentry, linewidth =0.5, rstride=5, cstride=125, color = color)
ax.yaxis.set_major_formatter(LogFormatter())
ax.yaxis.set_major_locator(plt.MaxNLocator(maxYOrderOfTen))
ax.set_yticklabels(myLogAxesList)
pos = ax.get_position()
pos.x0 = pos.x0 - 0.01
pos.x1 = pos.x1 - 0.01
ax.set_position(pos)

ax.set_ylabel(r"Concentration ($\mu$M)")
ax.set_xlabel('Voltage (mV)')
ax.set_zlabel('Entry Rate (molecules/s)')

ax.tick_params(axis= 'both')
axarr[1,3].axis('off')


#=============================================================================
#=========================== Final print statements ==========================
#=============================================================================

print("\n")
print("Model results: \n")
print("The entry rate for", Nch, "channel(s) at 1 uM and", 
      voltage_For_Model, "V is {:3.4f}".format(Nentry1) + " molecules per second.")
print("The entry rate for", Nch, "channel(s) at 2 uM and", 
      voltage_For_Model, "V is {:3.4f}".format(Nentry2) + " molecules per second.")
print("The entry rate for", Nch, "channel(s) at 3 uM and", 
      voltage_For_Model, "V is {:3.4f}".format(Nentry3) + " molecules per second.")
print("The entry rate for", Nch, "channel(s) at 100 uM and", 
      voltage_For_Model, "V is {:3.4f}".format(Nentry100) + " molecules per second.")
print("The entry rate for", Nch, "channel(s) at 10000 uM and", 
      voltage_For_Model, "V is {:3.4f}".format(Nentry10000) + " molecules per second.")

print("The maximum entry rate is {}".format(round(Max_entry,3)), "at {}".format(round(Max_entry_concentration,3)), "uM")
print("The half-entry rate is is {}".format(Half_entry), "at {}".format(round(Half_entry_concentration, 3)), "uM")
print("The Potential of Max Block is {:3.4f}".format(V0 *1000)+ " mV.", "CI: {:3.4f}".format(V0_CI))
print("The E1 is {:3.4}".format(E1), "CI: {}".format(E1_CI))
print("The E2 is {:3.4}".format(E2), "CI: {:3.4f}".format(E2_CI))
print("The k_1 is {:3.4}".format(k_1), "CI: {:3.4f}".format(k_1_CI))
print("The k2 is {:3.4}".format(k2), "CI: {:3.4f}".format(k2_CI))
print("The k_2 is {:3.4}".format(k_2), "CI: {:3.4f}".format(k_2_CI))


# make it full-screen!
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

