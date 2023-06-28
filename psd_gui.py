# -*- coding: utf-8 -*
"""
Created on Mon May 13 09:15:18 2019 by Jakob Weyel, Eduard-Zintl-Institut für
Anorganische und Physikalische Chemie, TU Darmstadt

Make sure that your graphics backend is set to 'Tkinter' for functions such as
'PeakPicking' and 'Show_Peaks'.

If want, you can include a fourier series of your choice (change the parameter
k as you like) in 'PSD_calc' instead of a simple sin function to describe the
periodic stimulation which is part of the convolution in the fourier
transformation.

@author: jakob.weyel@tu-darmstadt.de
"""

import os
import time
import re

from tkinter import *
from tkinter import messagebox
from tkinter import filedialog

import numpy as np
from scipy import integrate as igr
import matplotlib.pyplot as plt
import pandas as pd

# The following two lines import stylesheets to format graphs. If you don't have any, comment them out
#if os.environ['LOGNAME'] == 'jakub' :
#    plt.style.use('/home/jakub/HESSENBOX-DA/Diverses/TU_Design.mplstyle')

'''
_______________________________________________________________________________
Dialogue box (YES/NO) with save button
_______________________________________________________________________________
'''

def yesno(name, output, text): # Decision if something will be saved
    
    msgbox = messagebox.askquestion ('Obacht!', text, icon = 'warning')
    if msgbox == 'yes':
       output.to_csv(name, sep = '\t', index = False)
       messagebox.showinfo('Yeah, man!', 'Saved as: ' + name)
    else:
        messagebox.showinfo('Allright.','Nothing got saved.')

'''
_______________________________________________________________________________
Open files
_______________________________________________________________________________
'''

def FileOpen(text): # Loads the desired data
	filename = filedialog.askopenfilename(title = text)
	return filename

'''
_______________________________________________________________________________
functions
_______________________________________________________________________________
'''

def PSD_calc(): # Calculates PSD spectra
    # Number of spectra per period, periods to cut off, phase resolution and the path are input via the GUI
    
    start = time.time()
    
    n_sp = int(Entry_n_sp.get()) # Number of spectra per period
    cutoff_per = int(Entry_cutoff_per.get()) #Number of periods to cut off
    cutoff_sp = n_sp*cutoff_per # Calculated number of spectra to cut off
    dphi = int(Entry_dphi.get()) # Phasse resolution
    
    # Reading data of different instruments/software having different formats    
    if instrument.get() == "Bruker/OPUS (DRIFTS)":
        text = 'Path of your catalyst spectra'
        name_dataCat = FileOpen(text)
        
        name_t = name_dataCat.split('.')
        name_t = name_t[0] + '_t.' + name_t[1]
        
        text = 'Path of your reference spectra'    
        name_dataRef = FileOpen(text)
        
        # Data of catalyst spectra and reference spectra
        dataCat = pd.read_csv(r''+name_dataCat, sep="\t", header = None)
        dataCat = dataCat.values
        t_inp = np.genfromtxt(r''+name_t, delimiter="\t")
        
        if name_dataRef != '':
            dataRef = pd.read_csv(r''+name_dataRef, sep="\t", header = None)
            dataRef = dataRef.values
            
            # Calculating the difference for gas phase subtraction or whatever you want
            data = dataCat-dataRef
            data[:,0] = dataCat[:,0]
            
        else :
            data = dataCat
        
    elif instrument.get() == "Horiba/LabSpec (Raman)":
        text = 'Path of your catalyst spectra'
        name_dataCat = FileOpen(text)
        
        text = 'Path of your reference spectra'
        name_dataRef = FileOpen(text)
        
        # Data of catalyst spectra and reference spectra
        dataCat = pd.read_csv(r''+name_dataCat, sep="\t", header = None)
        dataCat = dataCat.values
        
        # Calculating t_inp into seconds
        t_inp = dataCat[1:,0]
        t_inp = np.reshape(t_inp,(t_inp.size,1))
        t_inp = (t_inp - t_inp[0]) * 24 * 3600
        t_inp = t_inp + t_inp[1]
        
        dataCat = np.delete(dataCat, 0, axis=1)
        dataCat = dataCat.T
        
        if name_dataRef != '':
            dataRef = pd.read_csv(r''+name_dataRef, sep="\t", header = None)
            dataRef = dataRef.values
            dataRef = np.delete(dataRef, 0, axis=1)
            dataRef = dataRef.T
            
            # Calculating the difference for gas phase subtraction or whatever you want
            data = dataCat-dataRef
            data[:,0] = dataCat[:,0]
            
        else :
            data = dataCat
                  
    # If t_inp is 1D array convert to 2D array for proper indexing
    #if t_inp.ndim == 1:
    #    t_inp = np.reshape(t_inp,(t_inp.size,1))
    
       
    if cutoff_sp != 0:
        # Cut off spectra from the beginning
        data = np.delete(data, np.s_[1:cutoff_sp+1], axis = 1)
        t_inp = np.delete(t_inp,np.s_[-(cutoff_sp):],axis = 0)
        
    
    n_per = t_inp.shape[0]/n_sp # number of periods by dividing number of spectra by number of spectra per period
    
    
    # Averaging all periods into one period
    
    Energy_values = (data[:,0]) # Cache the energy values / wavenumbers
    
    Energy_values = np.reshape(Energy_values,(Energy_values.size,1)) # Make 2D array for further computations
    
    spectra_per = np.split(data[:,1:], n_per, axis = 1) # Split the wholeness of all spectra in 'data' into minor ndarrays for each period
    
    sum_spectra_per = np.divide(sum(spectra_per),n_per) # sum up all cells of the created ndarrays that have the same index and divide by the number of periods
        
    data = np.concatenate((Energy_values, sum_spectra_per),axis = 1) # Concatenate Energy_values and averaged spectra back into 'data'
    
    # Defining the 0 of the intensity axis as the mean value of all time resolved datapoints for each particular wavenumber
    
    i = 0
    for row in data:
        data[i,1:] = data[i,1:]-np.mean(data[i,1:])
        i = i+1
        
    # Definition of values for the fourier transformation
    t_per = t_inp[n_sp-1,0]
    omega = 2*np.pi/t_per # omega in s^-1
    phi = np.arange(0,361,dphi) # phi is the phase shift which occurs as answer of the system to the external stimulation in the experiment
    
    spectra = np.zeros((len(data[:,0]),len(phi)+1))
    spectra[:,0] = dataCat[:,0]
    dummy = spectra
    
    # Do the fourier transformation for all predefined values of phi
    for k in np.arange(1,2): # set k>1 for modeling a rectangular function via fourier synthesis which will be folded with the time resolved spectra
        for i in range(1,len(phi)+1):
            for j in range(0,len(data[:,0])):
                dummy[j,i] = 2/t_inp[int(n_sp),0]*igr.trapz(data[j,1:]*(1/(2*k))*np.sin((2*k-1)*omega*t_inp[0:n_sp,0]+phi[i-1]*2*np.pi/360))
                
    spectra[:,1:] = spectra[:,1:]+dummy[:,1:]
    
    # Plot spectra (if needed, uncomment it)
    # plt.figure()
    # for i in range(1,len(phi)+1):
    
    #     plt.plot(spectra[:,0],spectra[:,1:],label = str(phi[i-1]))
    
    # plt.ylim(-np.amax(spectra[:,1:]), np.amax(spectra[:,1:]))
    
    # plt.xlabel(xUnit_list[xUnit.get()]) #Gets x axis label
    # plt.ylabel(yUnit_list[yUnit.get()]) #Gets y axis label
    
    # add legend
    # plt.legend(phi, title = r'$\varphi$ / °', loc = 'upper right')
        
    # Put PSD spectra into data frame
    phi_str = [str(item) + '°' for item in phi]
        
    header = ['Wavenumber'] + phi_str
    output = pd.DataFrame(data = spectra[:,:], columns = header)
    
    end = time.time()
    
    # Save the PSD spectra
    
    text = 'Shall the PSD spectra be saved as .txt?'
    name = name_dataCat.split('.')
    name = name[0] + '_PSD_spectra_' + str(cutoff_per) + '_periods_cutoff_' + str(dphi) + '_dphi.txt'
    yesno(name, output, text)
    print('Runtime PSD: ' + str(end-start) +' s.')

    return

def Spectra_diff():
    # Loads number of spectra per period from textbox 'Entry_n_sp'
    text = 'Path of your catalyst spectra'
    name_dataCat = FileOpen(text)
    
    text = 'Path of your reference spectra'
    name_dataRef = FileOpen(text)
    
    # Reading data of different instruments/software having different formats
    if instrument.get() == "Bruker/OPUS (DRIFTS)":
        # Data of catalyst spectra and reference spectra
        dataCat = pd.read_csv(r''+name_dataCat, sep="\t", header = None)
        dataCat = dataCat.values
        
        dataRef = pd.read_csv(r''+name_dataRef, sep="\t", header = None)
        dataRef = dataRef.values
            
    elif instrument.get() == "Horiba/LabSpec (Raman)":
        # Data of catalyst spectra and reference spectra
        dataCat = pd.read_csv(r''+name_dataCat, sep="\t", header = None)
        dataCat = dataCat.values
        dataCat = np.delete(dataCat, 0, axis=1)
        dataCat = dataCat.T
        
        dataRef = pd.read_csv(r''+name_dataRef, sep="\t")#, header = None)
        dataRef = dataRef.values
        dataRef = np.delete(dataRef, 0, axis=1)
        dataRef = dataRef.T
        
    # Data of catalyst spectra and reference spectra
    dataCat = pd.read_csv(r''+name_dataCat, sep="\t", header = None) # if dataset got header: comment last argument out (can't compute strings!))
    dataCat = dataCat.values
    
    dataRef = pd.read_csv(r''+name_dataRef, sep="\t") # , header = None)
    dataRef = dataRef.values
    
    # Calculating the difference for gas phase subtraction or whatever you want
    data = dataCat-dataRef
    data[:,0] = dataCat[:,0]
    
    # plots the difference spectra
    # plt.figure()
    
    # plt.plot(data[:,0],data[:,1:])
    
    # # grabs x and y units for graph
    # plt.xlabel(xUnit_list[xUnit.get()]) # Gets x axis label
    # plt.ylabel(yUnit_list[yUnit.get()]) # Gets y axis label
    
    # plt.ylim(-np.amax(data[:,1:]), np.amax(data[:,1:]))
    
    # Put difference spectra into data frame
    diff = ['diff ' + str(i) for i in range(1,len(data[0,:]))]
    
    header = ['Wavenumber'] + diff
    output = pd.DataFrame(data = data[:,:], columns = header)
    
    # Save data frame to a file?
    
    text = 'Shall the difference spectra be saved as .txt?'
    name = name_dataCat.split('.') #Will be used as filename
    name = name[0] + '_diff_spectra.txt'
    yesno(name, output, text)

    return

def PointPicking():
    fig = plt.figure()
    text = 'Path of your spectra'
    name_dataSpectra = FileOpen(text)

    if 'PSD' in name_dataSpectra:
        # Data of PSD spectra
        dataSpectra = pd.read_csv(r''+name_dataSpectra, sep="\t")
        dataSpectra = dataSpectra.values
        # Plot all spectra at once
        plt.clf()
        plt.xlabel(xUnit_list[xUnit.get()]) # Gets x axis label
        plt.ylabel(yUnit_list[yUnit.get()]) # Gets y axis label
        
        for i in range(1, dataSpectra.shape[1]):
            plt.plot(dataSpectra[0:,0],dataSpectra[0:,i], 'o', picker = 3)
            plt.ylim(-np.amax(dataSpectra[:,1:]), np.amax(dataSpectra[:,1:]))
            plt.show()
            
            bands = []
            
            def onpick(event): # This function allows you to save the current position of your mouse to an array when clicked on a point of the shown graph
                thisline = event.artist
                xdata = thisline.get_xdata()
                ind = event.ind
                points = xdata[ind]
                bands.append(points[int(len(points)/2)]) # Append currently clicked on position to array
                            
                # Save as soon as a point is clicked on
                bands_new = sorted(set(bands)) # sorted sorts, set sorts and removes doubly counted ones
                name = name_dataSpectra.split('.')
                bands_new = np.array(bands_new)
                np.savetxt(name[0] + '_peaks.txt',bands_new, delimiter = '\t')
                
        fig.canvas.mpl_connect('pick_event', onpick)
        
    else:
        n_sp = int(Entry_n_sp.get()) # Number of spectra per period
        cutoff_per = int(Entry_cutoff_per.get()) #Number of periods to cut off
        cutoff_sp = n_sp*cutoff_per # Calculated number of spectra to cut off
        
        # Reading data of different instruments/software having different formats
        if instrument.get() == "Bruker/OPUS (DRIFTS)":
            # Data of catalyst spectra
            dataSpectra = pd.read_csv(r''+name_dataSpectra, sep="\t")
            dataSpectra = dataSpectra.values
            
            name_t = name_dataSpectra.split('.')
            name_t = name_t[0] + '_t.' + name_t[1]        
            t_inp = np.genfromtxt(r''+name_t, delimiter="\t")
                
        elif instrument.get() == "Horiba/LabSpec (Raman)":
            # Data of catalyst spectra
            dataSpectra = pd.read_csv(r''+name_dataSpectra, sep="\t", header = None)
            dataSpectra = dataSpectra.values
            
            # Calculating t_inp into seconds
            t_inp = dataSpectra[1:,0]
            t_inp = np.reshape(t_inp,(t_inp.size,1))
            t_inp = (t_inp - t_inp[0]) * 24 * 3600
            t_inp = t_inp + t_inp[1]
            
            dataSpectra = np.delete(dataSpectra, 0, axis=1)
            dataSpectra = dataSpectra.T
       
        if cutoff_sp != 0:
            # Cut off spectra from the beginning
            dataSpectra = np.delete(dataSpectra, np.s_[1:cutoff_sp+1], axis = 1)
            t_inp = np.delete(t_inp,np.s_[-(cutoff_sp):],axis = 0)
            
        n_per = t_inp.shape[0]/n_sp # number of periods by dividing number of spectra by number of spectra per period
        
        # Averaging all periods into one period        
        Energy_values = (dataSpectra[:,0]) # Cache the energy values / wavenumbers
        
        Energy_values = np.reshape(Energy_values,(Energy_values.size,1)) # Make 2D array for further computations
        
        spectra_per = np.split(dataSpectra[:,1:], n_per, axis = 1) # Split the wholeness of all spectra in 'data' into minor ndarrays for each period
        
        sum_spectra_per = np.divide(sum(spectra_per),n_per) # sum up all cells of the created ndarrays that have the same index and divide by the number of periods
            
        dataSpectra = np.concatenate((Energy_values, sum_spectra_per),axis = 1) # Concatenate Energy_values and averaged spectra back into 'data'

        # Plot a view spectra at once
        plt.clf()
        plt.xlabel(xUnit_list[xUnit.get()]) # Gets x axis label
        plt.ylabel(yUnit_list[yUnit.get()]) # Gets y axis label
        
        n_plot = 10 # number of spectra to plot
        j = n_sp/n_plot
        if divmod(n_sp,n_plot)[1] != 0:
            j = (n_sp - divmod(n_sp,n_plot)[1])/n_plot
        
        #for i in range(1, dataCat.shape[1]-1):
        for i in range(1, n_plot+1):
            plt.plot(dataSpectra[0:,0],dataSpectra[0:,int(i*j)], 'o', picker = 3)
            plt.ylim(-np.amax(dataSpectra[:,1:]), np.amax(dataSpectra[:,1:]))
            plt.show()
            
            bands = []
            
            def onpick(event): # This function allows you to save the current position of your mouse to an array when clicked on a point of the shown graph
                thisline = event.artist
                xdata = thisline.get_xdata()
                ind = event.ind
                points = xdata[ind]
                bands.append(points[int(len(points)/2)]) # Append currently clicked on position to array
                            
                # Save as soon as a point is clicked on
                bands_new = sorted(set(bands)) # sorted sorts, set sorts and removes doubly counted ones
                name = name_dataSpectra.split('.')
                bands_new = np.array(bands_new)
                np.savetxt(name[0] + '_BLpoints.txt',bands_new, delimiter = '\t')
                
        fig.canvas.mpl_connect('pick_event', onpick)
    
    return

def Baseline(): # generate and substract baseline and PSD
    n_sp = int(Entry_n_sp.get()) # Number of spectra per period
    cutoff_per = int(Entry_cutoff_per.get()) #Number of periods to cut off
    cutoff_sp = n_sp*cutoff_per # Calculated number of spectra to cut off
    
    text = 'Path of your catalyst spectra'
    name_dataCat = FileOpen(text)
    
    text = 'Path of your baseline points'
    name_points = FileOpen(text)
    
    name_t = name_dataCat.split('.')
    name_t = name_t[0] + '_t.' + name_t[1]
    
    # Data of catalyst spectra
    dataCat = pd.read_csv(r''+name_dataCat, sep="\t", header = None)
    dataCat = dataCat.values
    t_inp = np.genfromtxt(r''+name_t, delimiter="\t")
    
    if cutoff_sp != 0:
        # Cut off spectra from the beginning
        dataCat = np.delete(dataCat, np.s_[1:cutoff_sp+1], axis = 1)
        t_inp = np.delete(t_inp,np.s_[-(cutoff_sp):],axis = 0)
        
    
    n_per = t_inp.shape[0]/n_sp # number of periods by dividing number of spectra by number of spectra per period
    
    
    # Averaging all periods into one period
    
    Energy_values = (dataCat[:,0]) # Cache the energy values / wavenumbers
    
    Energy_values = np.reshape(Energy_values,(Energy_values.size,1)) # Make 2D array for further computations
    
    spectra = dataCat[:,1:]    
    
    # generating baseline
    point_pos = np.genfromtxt(r''+name_points, delimiter="\n")
    # peak_pos = peak_pos[::-1] #Needs to be inverted otherwise it is the wrong way around!
    
    '''Compares every value in peak_pos with the wavenumbers from psd_spectra
    and the closest value is taken'''
    index_Energy = np.zeros(len(point_pos))
    i = 0
    for val in point_pos:
        point_pos[i] = min(Energy_values, key=lambda x:abs(x-val))        
        index_Energy[i] = abs(Energy_values-point_pos[i]).argmin()
        #Energy_values.index(point_pos[i])   
        i = i+1
        
    index_Energy = np.flip(index_Energy)    

    n_avg = 5 # number of data points to average each intensity (must be odd)

    I_mean = np.zeros((int(n_sp*n_per), len(point_pos)))
    I_bl = np.zeros((len(spectra), int(n_sp*n_per)))
    m = np.zeros((int(n_sp*n_per), len(point_pos)+1))
    b = np.zeros((int(n_sp*n_per), len(point_pos)+1))
    
    for i in range(len(point_pos)+1):
        if i == 0:
            I_mean_start = np.mean(spectra[:n_avg,:], axis=0)
            I_mean[:,i] = np.mean(spectra[int(index_Energy[i]-(n_avg-1)/2):int(index_Energy[i]+(n_avg-1)/2+1), :], axis=0)
            m[:,i] = (I_mean[:,i]-I_mean_start)/(Energy_values[int(index_Energy[i])]-Energy_values[0])
            b[:,i] = I_mean_start - m[:,i] * Energy_values[0]
            I_bl[:int(index_Energy[i]),:] = m[:,i] * Energy_values[:int(index_Energy[i])] + b[:,i]
        elif i > 0 and i < len(point_pos):
            I_mean[:,i] = np.mean(spectra[int(index_Energy[i]-(n_avg-1)/2):int(index_Energy[i]+(n_avg-1)/2+1), :], axis=0)
            m[:,i] = (I_mean[:,i]-I_mean[:,i-1])/(Energy_values[int(index_Energy[i])]-Energy_values[int(index_Energy[i-1])])
            b[:,i] = I_mean[:,i] - m[:,i] * Energy_values[int(index_Energy[i])]
            I_bl[int(index_Energy[i-1]):int(index_Energy[i]),:] = m[:,i] * Energy_values[int(index_Energy[i-1]):int(index_Energy[i])] + b[:,i]
        else:
            I_mean_end = np.mean(spectra[len(Energy_values)-n_avg:len(Energy_values),:], axis=0)
            m[:,i] = (I_mean_end-I_mean[:,i-1])/(Energy_values[len(Energy_values)-1]-Energy_values[int(index_Energy[i-1])])
            b[:,i] = I_mean_end - m[:,i] * Energy_values[len(Energy_values)-1]
            I_bl[int(index_Energy[i-1]):len(Energy_values),:] = m[:,i] * Energy_values[int(index_Energy[i-1]):len(Energy_values)] + b[:,i]
   
    data_baseline = np.concatenate((Energy_values, I_bl),axis = 1) # Concatenate Energy_values and baseline back into 'data'
    output = pd.DataFrame(data = data_baseline[:,:], columns = None)    
    
    # Save the baseline
    
    text = 'Shall the baseline be saved as .txt?'
    name = name_dataCat.split('.')
    name = name[0] + '_baseline.txt'
    yesno(name, output, text)    
    
    return

def Show_Peaks():
    text = 'Path of your peaks'
    name_peaks = FileOpen(text)
    peaks = np.genfromtxt(r''+name_peaks, delimiter="\n")
    for ii in np.arange(0,len(peaks)):
        ymin,ymax = plt.gca().get_ylim()
        plt.plot([peaks[ii],peaks[ii]],[ymin,ymax],'r', linewidth = 0.5)
    
    return

def in_phase_angle():
    # Input: PSD spectra and peak positions
    text = 'Path of your PSD spectra'
    name_psd = FileOpen(text)    
    psd_spectra = pd.read_csv(r''+name_psd, delimiter="\t")
    # Needs to round because pandas uses more decimals than numpy
    psd_spectra.Wavenumber = pd.Series([round(val, 5) for val in psd_spectra.Wavenumber], index = psd_spectra.index)
    
    text = 'Path of your peaks'
    name_peaks = FileOpen(text)
    peak_pos = np.genfromtxt(r''+name_peaks, delimiter="\n")
    # peak_pos = peak_pos[::-1] #Needs to be inverted otherwise it is the wrong way around!
    
    '''Compares every value in peak_pos with the wavenumbers from psd_spectra
    and the closest value is taken'''
    i = 0
    for val in peak_pos:
        peak_pos[i] = min(psd_spectra.Wavenumber, key=lambda x:abs(x-val))
        i = i+1
    
    # Reading time values of different instruments/software having different formats
    if instrument.get() == "Bruker/OPUS (DRIFTS)":
        text = 'Path of your time values'
        name_t = FileOpen(text)
        t_inp = np.genfromtxt(r''+name_t, delimiter="\t")
            
    elif instrument.get() == "Horiba/LabSpec (Raman)":
        text = 'Path of your catalyst spectra because of needed time values'
        name_dataCat = FileOpen(text)
        dataCat = pd.read_csv(r''+name_dataCat, sep="\t", header = None)
        dataCat = dataCat.values     
        
        # Calculating t_inp into seconds
        t_inp = dataCat[1:,0]
        t_inp = np.reshape(t_inp,(t_inp.size,1))
        t_inp = (t_inp - t_inp[0]) * 24 * 3600
        t_inp = t_inp + t_inp[1]
    
    # If t_inp is 1D array convert to 2D array for proper indexing
    #if t_inp.ndim == 1:
    #    t_inp = np.reshape(t_inp,(t_inp.size,1))
        
    n_sp = int(Entry_n_sp.get())
    t_per = t_inp[n_sp-1,0]
    
    # Separate the rows belonging to the chosen peak positions 
    phi_at_peaks = psd_spectra[psd_spectra.Wavenumber.isin(peak_pos)]
    phi_at_peaks = phi_at_peaks.iloc[:,1:] #Get rid of the wavenumber column
    
    # Read out the phase angle belonging to the respective maximum
    Wmax = phi_at_peaks.idxmax(axis = 1)
    Wmax = np.array(Wmax.values,dtype = str)

    numbers = re.compile(r'\d+(?:\.\d+)?') # define that only characters important for decimals are kept
    Wmax = np.array( [numbers.findall(item) for item in Wmax] ).reshape(-1) # finds only numbers and dots and puts them into the array

    # Wmax = np.core.defchararray.replace(Wmax,'°','')
    Wmax = np.array(Wmax,dtype = int)
    
    # Convert phase angle at maximum into time at maximum
    tmax = (360-Wmax)/360*t_per
    
    # Rounds wavenumber and time value to one decimal and puts all into data frame
    peak_pos = np.around(peak_pos,1)
    tmax = np.around(tmax,1)
    
    output = pd.DataFrame({'Wavenumber': peak_pos, 'Phi_max / °': Wmax, 't / s with t_Per. = ' + str(t_per) + ' s': tmax})
    
    # Save dataframe
    text = 'Shall peak positions and in phase angles be saved as .txt?'
    name = name_psd.split('.') #will be used as filename in case of saving
    name = name[0] + '_peaks_iPW.txt'
    yesno(name, output, text)
  
def Show_Graph(): # Plots any graph you want
    text = 'Path of your spectra'
    name_psd = FileOpen(text)
    psd = pd.read_csv(r''+name_psd, sep="\t")
    psd = psd.values
    
    # Plot all spectra at once
    plt.figure()
    
    
    # plt.plot(psd[:,0],psd[:,121:142:2]) # takes one period out of TR spectra and plots only 11 of them (more looks crowded)
    # plt.plot(psd[:,0],psd[:,241:322:8]) # takes one period out of TR spectra and plots only 11 of them (more looks crowded)
    # plt.plot(psd[:,0],psd[:,1201:1442:22]) # takes one period out of TR spectra and plots only 11 of them (more looks crowded)
    plt.plot(psd[:,0],psd[:,1:])
    
    # grabs x and y units for graph 
    plt.xlabel(xUnit_list[xUnit.get()]) # Gets x axis label
    plt.ylabel(yUnit_list[yUnit.get()]) # Gets y axis label
    
    # plt.yticks([],[])
    plt.ylim(np.amin(psd[:,1:]), np.amax(psd[:,1:]))
    plt.xlim(np.amin(psd[:,0]), np.amax(psd[:,0]))
    
    # phi = np.arange(0,360,30)
    # plt.legend(phi, title = r'$\varphi$ / °', loc = 'upper right')
    
    return

def contour():
    
    text = 'Path of your PSD spectra'
    name_psd = FileOpen(text)
    psd_spectra = pd.read_csv(r''+name_psd, sep='\t')
    psd_spectra = psd_spectra.values
    
    text = 'Path of your time values'
    name_t = FileOpen(text)
    t_inp = np.genfromtxt(r''+name_t, delimiter="\t")
    
    # If t_inp is 1D array convert to 2D array for proper indexing
    if t_inp.ndim == 1:
        t_inp = np.reshape(t_inp,(t_inp.size,1))
    
    wavenumber = psd_spectra[1:,0] # all wavenumbers
    
    spec = psd_spectra[1:,1:] # all intensities

    spec = spec[:,::-1] # inverse for correct spectrum-time relation
    
    n_sp = int(Entry_n_sp.get()) # number of spectra per period
    
    t_per = t_inp[n_sp-1,0] # period length
    
    t = np.linspace(0,t_per, spec.shape[1]) # time vector
    
    WN, T = np.meshgrid(wavenumber, t) # create 2D grid
    
    # plot stuff
    
    plt.figure()
    
    plt.contourf(WN, T, -spec.T, 100, cmap = 'gist_heat')
    
    plt.xlabel(xUnit_list[xUnit.get()])
    plt.ylabel(r'$t_\mathrm{\varphi}$ / s')
    
    return

def course():
    # Input: time resolved specta, time values, peak positions
    n_sp = int(Entry_n_sp.get()) #Number of spectra per period
    cutoff_per = int(Entry_cutoff_per.get()) #Number of periods to cut off
    cutoff_sp = n_sp*cutoff_per #Calculated number of spectra to cut off
    
    # text = 'Path of your catalyst spectra'
    # name_dataCat = FileOpen(text)
    
    # name_t = name_dataCat.split('.')
    # name_t = name_t[0] + '_t.' + name_t[1]
    
    # text = 'Path of your reference spectra'
    # name_dataRef = FileOpen(text)
    
    # text = 'Path of your peak positions'
    # name_peaks = FileOpen(text)        
    
    # Data of catalyst spectra, reference spectra and peak positions
    #dataCat = pd.read_csv(r''+name_dataCat, sep="\t", header = None)
    #dataCat = dataCat.values
    #t_inp = np.genfromtxt(r''+name_t, delimiter="\t")
    #peaks = np.genfromtxt(r''+name_peaks, delimiter="\n")
    
    #####
    # Reading data of different instruments/software having different formats    
    if instrument.get() == "Bruker/OPUS (DRIFTS)":
        text = 'Path of your catalyst spectra'
        name_dataCat = FileOpen(text)
        
        name_t = name_dataCat.split('.')
        name_t = name_t[0] + '_t.' + name_t[1]
        
        text = 'Path of your reference spectra'    
        name_dataRef = FileOpen(text)
        
        text = 'Path of your peak positions'
        name_peaks = FileOpen(text)
        peaks = np.genfromtxt(r''+name_peaks, delimiter="\n")
        
        # Data of catalyst spectra and reference spectra
        dataCat = pd.read_csv(r''+name_dataCat, sep="\t", header = None)
        dataCat = dataCat.values
        t_inp = np.genfromtxt(r''+name_t, delimiter="\t")
        
        if name_dataRef != '':
            dataRef = pd.read_csv(r''+name_dataRef, sep="\t", header = None)
            dataRef = dataRef.values
            
            # Calculating the difference for gas phase subtraction or whatever you want
            data = dataCat-dataRef
            data[:,0] = dataCat[:,0]
            
        else :
            data = dataCat
        
    elif instrument.get() == "Horiba/LabSpec (Raman)":
        text = 'Path of your catalyst spectra'
        name_dataCat = FileOpen(text)
        
        text = 'Path of your reference spectra'
        name_dataRef = FileOpen(text)
        
        text = 'Path of your peak positions'
        name_peaks = FileOpen(text)
        peaks = np.genfromtxt(r''+name_peaks, delimiter="\n")
        
        # Data of catalyst spectra and reference spectra
        dataCat = pd.read_csv(r''+name_dataCat, sep="\t", header = None)
        dataCat = dataCat.values
        
        # Calculating t_inp into seconds
        t_inp = dataCat[1:,0]
        t_inp = np.reshape(t_inp,(t_inp.size,1))
        t_inp = (t_inp - t_inp[0]) * 24 * 3600
        t_inp = t_inp + t_inp[1]
        
        dataCat = np.delete(dataCat, 0, axis=1)
        dataCat = dataCat.T
        
        if name_dataRef != '':
            dataRef = pd.read_csv(r''+name_dataRef, sep="\t", header = None)
            dataRef = dataRef.values
            dataRef = np.delete(dataRef, 0, axis=1)
            dataRef = dataRef.T
            
            # Calculating the difference for gas phase subtraction or whatever you want
            data = dataCat-dataRef
            data[:,0] = dataCat[:,0]
            
        else :
            data = dataCat
    #####
    
        
    # If t_inp is 1D array convert to 2D array for proper indexing
    #if t_inp.ndim == 1:
    #    t_inp = np.reshape(t_inp,(t_inp.size,1))
    
    
    if cutoff_sp != 0:
        # Cut off spectra from the beginning
        data = np.delete(data, np.s_[1:cutoff_sp+1], axis = 1)
        t_inp = np.delete(t_inp,np.s_[-(cutoff_sp):],axis = 0)
    
    t_1spectrum = t_inp[1,0]-t_inp[0,0] #Calculates the time needed for measuring one spectrum
    print('One spectrum needed ' + str(t_1spectrum) + ' s.')
    
    i = 0
    for val in peaks:
        peaks[i] = min(data[:,0], key=lambda x:abs(x-val)) #If the wavenumber from 'peaks' and data[:,0] don't fit 100% the value with the lowest deviation will be taken
        i = i+1
        
    # Plot the course
    pos = np.zeros(len(peaks))
    for i in np.arange(0,len(peaks)):
        if i%15 == 0: # opens a new plot window every XXX lines (insert number of your choice). Otherwise the colours get confusing
            # Highlight the different phases of your periodic stimulation
            plt.figure()
            if n_sp != 0:
                for j in np.arange(min(t_inp[:,0]),max(t_inp[:,0]/60),n_sp*t_1spectrum/60): #divide by 60 to get from s to min
                    plt.axvspan(j, j+n_sp/2*t_1spectrum/60, facecolor='k', alpha=0.25)
        
        dummy = np.where(data[:,0] == peaks[i])
        pos[i] = dummy[0]
        
        plt.plot(t_inp[:,0]/60,data[int(pos[i]),1:], label = str(np.around(peaks[i],0)))
        plt.legend(title = xUnit_list[xUnit.get()], loc='upper right')
        plt.xlabel('$t$ / min')
        plt.ylabel(yUnit_list[yUnit.get()]) # Gets y axis label
   
    # Write sth. to save output as txt
    output1 = pd.DataFrame({'t / min': t_inp[:,0]/60})
    
    output2=pd.DataFrame(data=data[pos.astype(int),1:].T, columns=np.around(peaks,0).astype(int))
    
    output = pd.merge(output1,output2, left_index=True, right_index=True)
    
    # Save data frame to a file?
    
    text = 'Shall the difference spectra be saved as .txt?'
    name = name_dataCat.split('.') #Will be used as filename
    name = name[0] + '_course.txt'
    yesno(name, output, text)
            
'''
_______________________________________________________________________________
GUI stuff
_______________________________________________________________________________
'''

PSD_GUI = Tk()
Entry_n_sp = StringVar()
Entry_cutoff_per = StringVar()
Entry_dphi = StringVar()

PSD_GUI.title('Have fun with PSD!')

Label_DD_instrument = Label(PSD_GUI, text = 'Choose your instrument/software!').pack()

instrument_list = ['Bruker/OPUS (DRIFTS)', 'Horiba/LabSpec (Raman)']
instrument = StringVar(PSD_GUI)
instrument.set('Bruker/OPUS (DRIFTS)')
DD_instrument = OptionMenu(PSD_GUI, instrument, *instrument_list)
DD_instrument.config(width=25)
DD_instrument.pack()

Label_n_sp = Label(PSD_GUI, text = 'Type in the number of spectra per period!').pack()

Entry_n_sp = Entry(PSD_GUI, textvariable = Entry_n_sp)
Entry_n_sp.insert(END,'320')
Entry_n_sp.pack()

Label_cutoff_per = Label(PSD_GUI, text = 'Choose the number of periods to cut off!').pack()

Entry_cutoff_per = Entry(PSD_GUI, textvariable = Entry_cutoff_per)
Entry_cutoff_per.insert(END,'0')
Entry_cutoff_per.pack()

Label_dphi = Label(PSD_GUI, text = 'Choose your phase resolution!').pack()

Entry_dphi = Entry(PSD_GUI, textvariable = Entry_dphi)
Entry_dphi.insert(END,'30')
Entry_dphi.pack()

Label_DD_yUnit = Label(PSD_GUI, text = 'Choose your y label!').pack()

yUnit_list = {'-lg(R)': r'-lg($R$)', 'Extinction': r'$E$', 'Transmission in %': r'$T$ / %', 'Absorption in %': r'$A$ / %', 'Intensity in a.U.': r'$I$ / a.U.'}
yUnit = StringVar(PSD_GUI)
if instrument.get() == "Bruker/OPUS (DRIFTS)":
    yUnit.set('-lg(R)')
elif instrument.get() == "Horiba/LabSpec (Raman)":
    yUnit.set('Intensity in a.U.')
DD_yUnit = OptionMenu(PSD_GUI, yUnit, *yUnit_list)
DD_yUnit.config(width=14)
DD_yUnit.pack()

Label_DD_xUnit = Label(PSD_GUI, text = 'Choose your x label!').pack()

xUnit_list = {'Wavenumber in 1/cm': r'$\tilde{\nu}$ / cm$^{-1}$', 'Energy in eV': r'$E$ / eV', 'Wavelength in nm': r'$\lambda$ / nm', 'Raman Shift in 1/cm': r'$\Delta\tilde{\nu}$ / cm$^{-1}$'}
xUnit = StringVar(PSD_GUI)
xUnit.set('Wavenumber in 1/cm')
DD_xUnit = OptionMenu(PSD_GUI, xUnit, *xUnit_list)
DD_xUnit.config(width=22)
DD_xUnit.pack()


Bt_PSD_calc = Button(PSD_GUI, text = 'TRS -> PSD', command = PSD_calc).pack()

Label_PointPicking = Label(PSD_GUI, text = 'If no graph is shown, close all tkinter windows except your graph! \n If a data point is clicked on, all data points \n clicked on until now are written into a file:').pack()

Bt_PointPicking = Button(PSD_GUI, text = 'point picking', command = PointPicking).pack()

Label_Baseline_corr = Label(PSD_GUI, text = 'Generate a baseline with picked points:').pack()

Bt_Baseline_corr = Button(PSD_GUI, text = 'baseline', command = Baseline).pack()

Label_Spectra_diff = Label(PSD_GUI, text = 'Here you calculate difference spectra:').pack()

Bt_Spectra_diff = Button(PSD_GUI, text = 'difference spectra', command = Spectra_diff).pack()

Label_Show_Peaks = Label(PSD_GUI, text = 'here you can highlight chosen peak positions:').pack()

Bt_Show_Peaks = Button(PSD_GUI, text = 'show peaks', command = Show_Peaks).pack()

Label_in_phase_angle = Label(PSD_GUI, text = 'Calculate in phase angle and in phase time:').pack()

Bt_in_phase_angle = Button(PSD_GUI, text = 'in phase angle', command = in_phase_angle).pack()

Label_Show_Graph = Label(PSD_GUI, text = 'Plot graphs you like:').pack()

Bt_Show_Graph = Button(PSD_GUI, text = 'show graph', command = Show_Graph).pack()

Label_contour = Label(PSD_GUI, text = 'Create contour plot from PSD spectra:').pack()

Bt_contour = Button(PSD_GUI, text = 'contour plot', command = contour).pack()

Label_course = Label(PSD_GUI, text = 'Create course plots of chosen spectra at chosen peak positions:').pack()

Bt_course = Button(PSD_GUI, text = 'course plot', command = course).pack()



PSD_GUI.mainloop()