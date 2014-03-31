import os
from xspec import *
import subprocess
from numpy import *


# Determinate the Hardness Ratio for a spectrum s1
def det_hr(s1):
    s1.ignore("**-0.5")
    s1.ignore("2.0-**")
    s = s1.rate[0]*s1.exposure
    s1.ignore("**-2.0")
    s1.notice("2.0-8.0")
    h = s1.rate[0]*s1.exposure
    HR = (h-s)/(h+s)
    s1.notice("0.5-8.0")
    errHR = sqrt((2*s*sqrt(h))**2+(2*h*sqrt(s))**2)/((h+s)**2)
    return HR, errHR


# Get the column density of the Milky Way
def getnh(ra, dec):
    output = subprocess.check_output("nh equinox=2000 ra="+str(ra)+" dec="
                                     + str(dec)+" | grep Weighted",
                                     shell=True)
    nhgal = float(output.split()[6])
    nhxspec = nhgal/10.**22
    return nhgal, nhxspec


# Plot the spectrum of the central source, with fit (zphabs*powerlw) and 
# residuals and gives total counts, rate, hardness ratiom, column density,
# photon index and reduced chi squared
def plot_spec(ID, name, ra, dec, z):

    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID[0:5]+'/repro')

    #Import spectrum
    s = Spectrum(ID+'_spec_grp.pi')
    counts = int(s.rate[0]*s.exposure)
    rate = s.rate[0]
    HR, errHR = det_hr(s)

    #Set model
    m = Model('phabs*zphabs*zpow')

    #Plotting setting
    Plot.xAxis = "keV"
    AllData.ignore('**-0.5, 8.0-**')
    if counts < 250:
        Plot.setRebin(0.1, 30)
    else:
        Plot.setRebin(10, 30)

    #nh
    nhgal, nhxspec = getnh(ra, dec)

    #Parameters
    m.phabs.nH = nhxspec
    m.phabs.nH.frozen = True
    m.zphabs.Redshift = z
    m.zpowerlw.Redshift = z
    m.zpowerlw.PhoIndex.values = [2, 0.01, 1.5, 1.7, 2.1, 2.3]

    #Fit
    Fit.statMethod = "cstat 1"
    Fit.query = "yes"
    Fit.nIterations = 50
    Fit.perform()

    #Strings
    nh = m.zphabs.nH.values[0]+m.phabs.nH.values[0]
    nherr = m.zphabs.nH.sigma+m.phabs.nH.values[0]
    photon = m.zpowerlw.PhoIndex.values[0]
    photonerr = m.zpowerlw.PhoIndex.sigma
    NH = 'NH = ('+str(nh)[0:4]+')e22'+'cm^-2'
    phot = 'Photon index = '+str(photon)[0:3]
    if Fit.dof == 0:
        redchi = 0
    else:
        redchi = Fit.statistic/Fit.dof
    hr = 'HR = '+str(HR)[0:4]
    chi = 'Red Chi^2 = '+str(redchi)[0:3]
    title = 'label top \"Observation ID: '+ID+' - '+name
    stuff = 'label f \"'+NH+'\n'+phot+'\n'+hr+'\n'+chi

    #Plot
    Plot.device = '/Users/andyscanzio/Documents/ETH/Semestri/'
    +'FS2014/Project/Plots/'+ID+'_counts.ps/cps'
    Plot.addCommand(title)
    Plot.addCommand(stuff)
    Plot('ldata residuals')
    Plot.device = 'none'

    AllData.clear()

    return counts, rate, HR, errHR, nh, nherr, photon, photonerr, redchi


# Plot the spectrum of the extended source, with fit (powerlw) and residuals
# and gives total counts, rate, photon index and reduced chi squared
def plot_spec_ext(ID, name, ra, dec, z):

    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID[0:5]+'/repro')

    #Import spectrum
    s = Spectrum(ID+'_spec_ext_grp.pi')
    counts = int(s.rate[0]*s.exposure)
    rate = s.rate[0]

    #Set model
    m = Model('phabs*zpow')

    #Plotting setting
    Plot.xAxis = "keV"
    AllData.ignore('**-0.5, 8.0-**')
    if counts < 250:
        Plot.setRebin(0.1, 30)
    else:
        Plot.setRebin(10, 30)

    #nh
    nhgal, nhxspec = getnh(ra, dec)

    #Parameters
    m.phabs.nH = nhxspec
    m.phabs.nH.frozen = True
    m.zpowerlw.Redshift = z

    #Fit
    Fit.statMethod = "cstat 1"
    Fit.query = "yes"
    Fit.nIterations = 50
    Fit.perform()

    #Strings
    photon = m.zpowerlw.PhoIndex.values[0]
    photonerr = m.zpowerlw.PhoIndex.sigma
    phot = 'Photon index = '+str(photon)[0:3]
    if Fit.dof == 0:
        redchi = 0
    else:
        redchi = Fit.statistic/Fit.dof
    chi = 'Red Chi^2 = '+str(redchi)[0:3]
    title = 'label top \"Observation ID: '+ID+' - '+name
    stuff = 'label f \"'+phot+'\n'+chi

    #Plot
    Plot.device = '/Users/andyscanzio/Documents/ETH/Semestri/'
    + 'FS2014/Project/Plots/'+ID+'_counts_ext.ps/cps'
    Plot.addCommand(title)
    Plot.addCommand(stuff)
    Plot('ldata residuals')
    Plot.device = 'none'

    AllData.clear()

    return counts, rate, photon, photonerr, redchi


# Write a txt file with the output of the plot_spec function
def output(ID, name, res):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/'
             + 'FS2014/Project/Outputs/')
    with open(ID+'_output.txt', 'w') as f:
        f.write('ObsID: '+ID+'\n')
        f.write('Name: '+name+'\n')
        f.write('Total count = '+str(res[0])+'\n')
        f.write('Counts per second = '+str(res[1])+'\n')
        f.write('Hardness ratio = '+str(res[2])[0:5]+' +/- '
                + str(res[3])[0:5]+'\n')
        f.write('NH = ('+str(res[4])[0:5]+' +/- '
                + str(res[5])[0:5]+')e22'+'cm^-2\n')
        f.write('Photon index = '+str(res[6])[0:5]+' +/- '
                + str(res[7])[0:5]+'\n')
        f.write('Reduced chi-squared = '+str(res[8]))
        f.write('\n')
        f.write('\n')


# Write a txt file with the output of the plot_spec_ext function
def output_ext(ID, name, res_ext):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/'
             + 'FS2014/Project/Outputs/')
    with open(ID+'_output_ext.txt', 'w') as f:
        f.write('ObsID: '+ID+'\n')
        f.write('Name: '+name+'\n')
        f.write('Total count = '+str(res_ext[0])+'\n')
        f.write('Counts per second = '+str(res_ext[1])+'\n')
        f.write('Photon index = '+str(res_ext[2])[0:5]+' +/- '
                + str(res_ext[3])[0:5]+'\n')
        f.write('Reduced chi-squared = '+str(res_ext[4]))
        f.write('\n')
        f.write('\n')


# Gives the hardness ration for an Observation
def HR(ID):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/'
             + 'FS2014/Project/Data/'+ID[0:5]+'/repro')
    s = Spectrum(ID+'_spec_grp.pi')
    HR, errHR = det_hr(s)
    AllData.clear()
    return HR, errHR


# Gives the rate of the extended source
def cps_ext(ID):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/'
             + 'FS2014/Project/Data/'+ID[0:5]+'/repro')
    s = Spectrum(ID+'_spec_ext_grp.pi')
    rate = s.rate[0]
    err_rate = 0
    return rate,  err_rate


# Gives the rate of the central source
def cps(ID):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/'
             + 'FS2014/Project/Data/'+ID[0:5]+'/repro')
    s = Spectrum(ID+'_spec_grp.pi')
    rate = s.rate[0]
    err_rate = 0
    return rate,  err_rate
