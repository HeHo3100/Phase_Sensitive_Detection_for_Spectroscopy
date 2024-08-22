# -*- coding: utf-8 -*
"""
Created on Mon May 13 09:15:18 2019 by Jakob Weyel, Eduard-Zintl-Institut für
Anorganische und Physikalische Chemie, TU Darmstadt

Edited since April 8 2024 by Henrik Hoyer, Eduard-Zintl-Institut für
Anorganische und Physikalische Chemie, TU Darmstadt

Make sure that your graphics backend is set to 'Tkinter' for functions such as
'PointPicking' and 'Show_Points'.

If want, you can include a fourier series of your choice (change the parameter
k as you like) in 'PSD_calc' instead of a simple sin function to describe the
periodic stimulation which is part of the convolution in the fourier
transformation.

@author: henrik.hoyer@tu-darmstadt.de
"""

import os
import time
import re

from tkinter import *
from tkinter import messagebox
from tkinter import filedialog
#from tkinter import Checkbutton
#from tkinter import IntVar

import numpy as np
from scipy import integrate as igr
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd


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

def PSD_calc():
    
    # input via the GUI 
    n_sp = int(Entry_n_sp.get()) # number of spectra per period
    n_per = int(Entry_n_per.get()) # number of periods that will be processed
    cutoff_per = int(Entry_cutoff_per.get()) # number of periods to cut off
    dphi = int(Entry_dphi.get()) # phase resolution [°]
    
    text = 'Path of your catalyst spectra' # choose data set
    name_data = FileOpen(text)
    
    start = time.time()
    
    
    # converts the DRIFTS data (OPUS, Bruker) into the right format
    if spectro.get() == "DRIFTS" and n_per != 1:
        
        # data set (wavenumber in first column; each column is one spectrum)
        data = pd.read_csv(r''+name_data, sep="\t", header = None)
        data = data.values
        
        # pick time values automatically 
        name_t = name_data.split('.')
        name_t = name_t[0] + '_t.' + name_t[1]
        
        t_inp = np.genfromtxt(r''+name_t, delimiter="\t") # time values [s]
        t_inp = np.delete(t_inp, np.s_[0,1,3,4], axis=1)
        
        energy_values = data[:,0] # wavenumbers [1/cm]    
        energy_values = np.reshape(energy_values,(energy_values.size,1)) # make 2D array for further computations
    
        data = np.delete(data, 0, axis=1) # delete energy values from data
        
    
    # converts the Raman data (LabSpec6, Horiba Jobin Yvon) into the right format
    elif spectro.get() == "Raman" and n_per != 1:
        
        # data set (time values in first column; Raman shift in first row; each row is one spectrum)
        data = pd.read_csv(r''+name_data, sep="\t", header=None)
        data = data.values
        
        # converts time values into the right unit
        t_inp = data[1:,0] # time values [0-1]
        t_inp = np.reshape(t_inp,(t_inp.size,1))
        t_inp = (t_inp - t_inp[0]) * 24 * 3600 # time values [s]
        
        t_dummy = t_inp[1,0] + t_inp[len(t_inp)-1,0] # total experiment time [s]
        t_dummy = t_dummy/(n_sp*n_per) # averaged aquisition time [s]
        for i in range(0,len(t_inp)):
            t_inp[i,0] = (i+1)*t_dummy

        energy_values = data[0,1:] # Raman shifts [1/cm]
        energy_values = np.reshape(energy_values,(energy_values.size,1))
        
        data = np.delete(data, 0, axis=1) # delete time values from data
        data = np.delete(data, 0, axis=0) # delete energy values from data
        data = data.T
        
        # remove cosmic ray spikes
        if Entry_spike_rem.get() == 1:
            spike_factor = 1.2 # spikes with x times higher intensity than in spectrum before or after are removed
            for i in range(1,len(data.T)-1): # from second to second last spectrum 
                for j in range(0,len(data)):
                    if data[j,i] > spike_factor*data[j,i-1] and data[j,i] > spike_factor*data[j,i+1]:
                        data[j,i] = (data[j,i-1]+data[j,i+1])/2

    
    # converts the UV-Vis data (AvaSoft, Avantes) into the right format   
    elif spectro.get() == "UV-Vis" and n_per != 1:
        
        # data set (time values in first row; wavelength in first column; each column is one spectrum)
        data = pd.read_csv(r''+name_data, sep=";", header=7)
        data = data.replace(',','.',regex=True)
        data = data.values

        # converts time values into the right unit
        t_inp = data[0,3:]  # time values [s]
        t_inp = np.reshape(t_inp,(t_inp.size,1))
        t_inp = t_inp.astype(float)
        t_inp[:,0] = t_inp[:,0] + t_inp[1,0]
        t_inp = t_inp/100000

        energy_values = data[1:,0] # wavelengths [nm]
        energy_values = np.reshape(energy_values,(energy_values.size,1))
        energy_values = energy_values.astype(float)

        data = np.delete(data, [0,1,2], axis=1) # delete energy values from data
        data = np.delete(data, 0, axis=0) # delete time values from data
        data = data.astype(float)
        
    
    # for the case that already averaged spectra are used
    elif n_per == 1:
        
        # data set (time values in first row; energy values in first column; each column is one spectrum)
        data = pd.read_csv(r''+name_data, sep="\t", header=0)
        data = data.values

        t_inp = data[0,1:] # time values [s]
        t_inp = np.reshape(t_inp,(t_inp.size,1))

        energy_values = data[1:,0] # Raman shift [1/cm]
        energy_values = np.reshape(energy_values,(energy_values.size,1))

        data = np.delete(data, 0, axis=1) # delete energy values from data
        data = np.delete(data, 0, axis=0) # delete time values from data
    
    
    # only for data that were not averaged before
    if n_per != 1:
        # cut off spectra at the end
        if n_per != len(data.T)/n_sp:
            cutoff_sp = (int(len(data.T)/n_sp) - n_per) * n_sp
            data = np.delete(data, np.s_[-(cutoff_sp):], axis=1)
            t_inp = np.delete(t_inp,np.s_[-(cutoff_sp):], axis=0)
        
        # cut off spectra from the beginning
        if cutoff_per != 0:
            cutoff_sp = cutoff_per * n_sp # Calculated number of spectra to cut off
            data = np.delete(data, np.s_[:cutoff_sp], axis=1)
            t_inp = np.delete(t_inp,np.s_[-(cutoff_sp):], axis=0)
          
        # average all periods into one period
        data = np.split(data, n_per-cutoff_per, axis=1) # split the wholeness of all spectra into minor ndarrays for each period    
        data = np.divide(sum(data),n_per) # sum up all cells of the created ndarrays that have the same index and divide by the number of periods       
        
        # bring energy_values in an ascending order and sort spectra correctly
        if energy_values[0,0] > energy_values[1,0]:
            energy_values = np.flip(energy_values)
            data = np.flip(data, axis=0)

    # normalize all averaged spectra to the single highest peak
    if spectro.get() == "Raman" and Entry_norm.get() == 1:
        I_max = np.max(data)
        # print("I_max = " + str(I_max)) # print the highest intensity to normalize other ranges
        data = data/I_max
        
    # only for data that were not averaged before (for saving)
    if n_per != 1:
        data_avg = np.concatenate((energy_values, data), axis=1) # concatenate energy values and averaged spectra
    
        t_avg = np.insert(t_inp[:n_sp,0], 0, 0) # insert 0 for later saving
        t_avg = np.reshape(t_avg,(t_avg.size,1))
        data_avg = np.concatenate((t_avg.T, data_avg), axis=0) # concatenate averaged spectra and time values
    
    # define the 0 of the intensity axis as the mean value of all time resolved datapoints for each particular energy value
    i = 0
    for row in data:
        data[i,:] = data[i,:]-np.mean(data[i,:])
        i = i+1
    
    # define the values for the phase sensitive detection
    t_per = t_inp[n_sp-1,0] # period length [s]
    omega = 2*np.pi/t_per # angular frequency of external stimulation [s]
    phi = np.arange(0,360,dphi) # phi is the phase shift which occurs as answer of the system to the external stimulation
    k_harmonic = int(Entry_k_harmonic.get()) # demodulation index to demodulate using a higher harmonic (1, 3, 5, ...) or a rectangular function (0)
    
    # phase sensitive detection for all predefined values of
    data_psd = np.zeros((len(data),len(phi)))
    for i in range(0,len(phi)):
        for j in range(0,len(data)): # if your external stimulation is more like a sine or a rectangular curve, comment the respective line out / in
            if k_harmonic == 0:
                data_psd[j,i] = 2/t_per*igr.trapz(data[j,:]*signal.square(omega*t_inp[:n_sp,0]+phi[i-1]*2*np.pi/360), x=t_inp[:n_sp,0]) # rectangular function
            else:
                data_psd[j,i] = 2/t_per*igr.trapz(data[j,:]*np.sin(k_harmonic*omega*t_inp[:n_sp,0]+phi[i-1]*2*np.pi/360), x=t_inp[:n_sp,0]) # sine curve
    
    data_psd = np.concatenate((energy_values, data_psd),axis = 1) # concatenate energy values and PSD spectra
    
    phi_conc = np.insert(phi, 0, 0) # insert 0 for later saving
    phi_conc = np.reshape(phi_conc,(phi_conc.size,1))
    data_psd = np.concatenate((phi_conc.T, data_psd), axis=0) # concatenate averaged spectra and phi values
                
    
    end = time.time()
    print('Runtime: ' + str(end-start) +' s.')
    
    
    # x and y unit of the spectra
    if spectro.get() == "DRIFTS":
        x_unit = 'Wavenumber / 1/cm'
        y_unit = '-log($R$) / a.u.'
    elif spectro.get() == "Raman":
        x_unit = 'Raman Shift / 1/cm'
        y_unit = 'Raman Intensity / a.u.'
    elif spectro.get() == "UV-Vis":
        x_unit = 'Wavelength / nm'
        y_unit = '-log($R$) / a.u.'
        
    # save the averaged spectra
    t_str = np.around(t_inp[:n_sp,0], 1)
    if n_per != 1:
        t_str = [' (' + str(item) + ' s)' for item in t_str]
        for i in range(0,n_sp):
            t_str[i] = '#' + str(i+1) + t_str[i]
        header = [x_unit] + t_str
        output = pd.DataFrame(data_avg, columns = header)
    
        text = 'Shall the averaged spectra be saved as .txt?'
        name = name_data.split('.')
        name = name[0] + '_' + str(cutoff_per) + '_periods_cutoff_' + 'averaged.txt'
        yesno(name, output, text)
    

    # save the PSD spectra
    phi_str = [str(item) + '°' for item in phi]       
    header = [x_unit] + phi_str
    output = pd.DataFrame(data_psd, columns = header)

    text = 'Shall the PSD spectra be saved as .txt?'
    name = name_data.split('.')
    name = name[0] + '_' + str(cutoff_per) + '_periods_cutoff_' + 'PSD_spectra_'  + str(dphi) + '_dphi.txt'
    yesno(name, output, text)
    
    
    
    # plot the PSD spectra (only if the number of spectra is not too big)
    if dphi >= 15:
        plt.figure(figsize=(10,5))
        plt.plot(data_psd[1:,0],data_psd[1:,1:])
        plt.xlabel(x_unit) # get x axis label
        plt.ylabel(y_unit) # get y axis label
        #plt.ylim(np.amin(data_psd[1:,1:]), np.amax(data_psd[1:,1:]))
        plt.xlim(np.amin(data_psd[1:,0]), np.amax(data_psd[1:,0]))
        plt.legend(phi, title = r'$\varphi$ / °', loc = 'upper left', bbox_to_anchor=(1,1)) # legend outside frame
    

    return

def Spectra_diff():
    
    text = 'Path of your averaged catalyst spectra' # choose data set
    name_data = FileOpen(text)
    
    text = 'Path of your averaged reference spectra' # choose reference data set (e.g., gas phase, baseline)
    name_dataRef = FileOpen(text)
    
    start = time.time()
    
    
    # data set (time values in first row; energy values in first column; each column is one spectrum)
    data = pd.read_csv(r''+name_data, sep="\t", header=0)
    data = data.values

    t_inp = data[0,1:] # time values [s]
    t_inp = np.reshape(t_inp,(t_inp.size,1))
    
    # reference dat set
    dataRef = pd.read_csv(r''+name_dataRef, sep="\t", header=0)
    dataRef = dataRef.values
    
    # calculate the difference of catalyst spectra and reference spectra
    data[1:,1:] = data[1:,1:] - dataRef[1:,1:]
    
    end = time.time()
    print('Runtime: ' + str(end-start) +' s.')
    
    
    # x and y unit of the spectra
    if spectro.get() == "DRIFTS":
        x_unit = 'Wavenumber / 1/cm'
        y_unit = '-log($R$) / a.u.'
    elif spectro.get() == "Raman":
        x_unit = 'Raman Shift / 1/cm'
        y_unit = 'Raman Intensity / a.u.'
    elif spectro.get() == "UV-Vis":
        x_unit = 'Wavelength / nm'
        y_unit = '-log($R$) / a.u.'
    
    t_str = np.around(t_inp[:,0], 1)
    t_str = [' (' + str(item) + ' s)' for item in t_str]
    for i in range(0,t_inp.size):
        t_str[i] = '#' + str(i+1) + t_str[i]
    header = [x_unit] + t_str
    output = pd.DataFrame(data, columns = header)
    
    text = 'Shall the difference spectra be saved as .txt?'
    name = name_data.split('.')
    name = name[0] + '_diff_spectra.txt'
    yesno(name, output, text)
    
    # plot the averaged spectra
    plt.figure(figsize=(10,5))
    plt.xlabel(x_unit) # get x axis label
    plt.ylabel(y_unit) # get y axis label
    #plt.ylim(np.amin(data[1:,1:]), np.amax(data[1:,1:]))
    plt.xlim(np.amin(data[1:,0]), np.amax(data[1:,0]))
    
    n_pl = 12 # maximum number of spectra to plot
    
    if n_pl < len(t_inp):
        t_legend = np.zeros((n_pl,1))
        i_pl = len(t_inp) // n_pl
        for i in range(n_pl):
            t_legend[i,0] = t_inp[(i+1)*i_pl-1,0]
            plt.plot(data[1:,0],data[1:,(i+1)*i_pl-1])
        plt.legend(t_legend[:n_pl,0], title = r'$t$ / s', loc = 'upper left', bbox_to_anchor=(1,1)) # legend outside frame
    
    else:
        plt.plot(data[1:,0],data[1:,1:])
        plt.legend(t_inp[:len(t_inp),0], title = r'$t$ / s', loc = 'upper left', bbox_to_anchor=(1,1)) # legend outside frame
    

    return

def PointPicking():
    
    text = 'Path of your spectra'
    name_data = FileOpen(text)
    data = pd.read_csv(r''+name_data, sep="\t", header=1)
    data = data.values
    
    # x and y unit of the spectra
    if spectro.get() == "DRIFTS":
        x_unit = 'Wavenumber / 1/cm'
        y_unit = '-log($R$) / a.u.'
    elif spectro.get() == "Raman":
        x_unit = 'Raman Shift / 1/cm'
        y_unit = 'Raman Intensity / a.u.'
    elif spectro.get() == "UV-Vis":
        x_unit = 'Wavelength / nm'
        y_unit = '-log($R$) / a.u.'


    # Plot all spectra at once
    fig = plt.figure(figsize=(10,5))
    plt.clf()
    plt.xlabel(x_unit) # Gets x axis label
    plt.ylabel(y_unit) # Gets y axis label
    
    
    n_pl = 12 # maximum number of spectra to plot
    
    if n_pl < data.shape[1]-1:
        i_pl = (data.shape[1]-1) // n_pl
        for i in range(n_pl):
            plt.plot(data[1:,0], data[1:,(i+1)*i_pl-1], 'o', picker=3)
            plt.show()
    else:
        plt.plot(data[1:,0],data[1:,1:], 'o', picker=3)
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
        name_points = name_data.split('.')
        bands_new = np.array(bands_new)
        np.savetxt(name_points[0] + '_points.txt',bands_new, delimiter = '\t')
            
    fig.canvas.mpl_connect('pick_event', onpick)
    
    return

def Show_Points():
    text = 'Path of your points'
    name_points = FileOpen(text)
    points = np.genfromtxt(r''+name_points, delimiter="\n")
    for ii in np.arange(0,len(points)):
        ymin,ymax = plt.gca().get_ylim()
        plt.plot([points[ii],points[ii]],[ymin,ymax],'r', linewidth = 0.5)
    
    return

def in_phase_angle():
    
    # x unit and of the spectra
    if spectro.get() == "DRIFTS":
        x_unit = 'Wavenumber / 1/cm'
    elif spectro.get() == "Raman":
        x_unit = 'Raman Shift / 1/cm'
    elif spectro.get() == "UV-Vis":
        x_unit = 'Wavelength / nm'
    
    # input: PSD spectra, point positions and averaged spectra for time values
    text = 'Path of your PSD spectra'
    name_psd = FileOpen(text)
    
    text = 'Path of your points'
    name_points = FileOpen(text)
    
    text = 'Path of your averaged spectra (for time values)'
    name_t = FileOpen(text)
    
    data_psd = pd.read_csv(r''+name_psd, delimiter="\t") # PSD spectra
    
    if name_points == '': # if no energy values are chosen, the in phase angles will be calculated for all energy values
        point_pos = data_psd.iloc[1:,0]
        point_pos = point_pos.values
        name = name_psd.split('.') # name for saving
        name = name[0] + '_iPW.txt'
    else:
        point_pos = np.genfromtxt(r''+name_points, delimiter="\n") # point positions
        point_pos = np.sort(point_pos) # sort points in ascending order
        name = name_points.split('.') # name for saving
        name = name[0] + '_iPW.txt'
    
    # compare every value in point_pos with the wavenumbers from psd_spectra and the closest value is taken
    i = 0
    for val in point_pos:
        point_pos[i] = min(data_psd[x_unit], key=lambda x:abs(x-val))
        i = i+1
    
    # read time values from averaged spectra data to convert the maximum phase angle into a time value
    t_inp = pd.read_csv(r''+name_t, delimiter="\t", header=0)
    t_inp = t_inp.values
    t_per = t_inp[0,len(t_inp.T)-1] # period length [s]
    
    # separate the rows belonging to the chosen point positions 
    phi_at_points = data_psd[data_psd[x_unit].isin(point_pos)]
    phi_at_points = phi_at_points.iloc[:,1:] # delete energy value in first cell
    
    # read out the phase angle belonging to the respective maximum
    w_max = phi_at_points.idxmax(axis = 1)
    w_max = np.array(w_max.values,dtype = str)

    numbers = re.compile(r'\d+(?:\.\d+)?') # define that only characters important for decimals are kept
    w_max = np.array( [numbers.findall(item) for item in w_max] ).reshape(-1) # find only numbers and dots and puts them into the array
    w_max = np.array(w_max,dtype = int)
    
    t_max = (360-w_max)/360*t_per # convert phase angle at maximum into time at maximum
    
    # generate negative t_max if intensity reaches maximum in second half of period
    t_max_neg = t_max - t_max + t_max # like this, t_max and t_max_neg are not connected in further computation
    for i in range(len(t_max_neg)):
        if t_max[i] > t_per/2:
            t_max_neg[i] = -(t_max[i] - t_per/2)
    
    # round wavenumbers and time values to one decimal and put all into data frame
    point_pos = np.around(point_pos,1)
    t_max_round = np.around(t_max,1)
    t_max_neg_round = np.around(t_max_neg,1)
    #output = pd.DataFrame({x_unit: point_pos, 'phi_max / °': w_max, 't / s with t_per = '+str(t_per)+' s': t_max_round})
    output = pd.DataFrame({x_unit: point_pos, 'phi_max / °': w_max, 't / s with t_per = '+str(t_per)+' s': t_max_round, '(-) t / s': t_max_neg_round})
    # save in phase angles and t_max
    text = 'Shall the in phase angles be saved as .txt?'
    # name = name_psd.split('.')
    # name = name[0] + '_points_iPW.txt'
    yesno(name, output, text)
    

    # plot the time constants
    plt.figure(figsize=(10,5))
    if len(point_pos) > 100: # if more than xxx time constants plot as connected lines
        plt.plot(point_pos, t_max_neg)
    else: # else plot as dots
        plt.plot(point_pos, t_max_neg, 'o')
    plt.xlabel(x_unit) # get x axis label
    plt.ylabel('$t$ / s') # get y axis label

    
    return
    
def Show_Graph():
    
    # x and y unit of the spectra
    if spectro.get() == "DRIFTS":
        x_unit = 'Wavenumber / 1/cm'
        y_unit = '-log($R$) / a.u.'
    elif spectro.get() == "Raman":
        x_unit = 'Raman Shift / 1/cm'
        y_unit = 'Raman Intensity / a.u.'
    elif spectro.get() == "UV-Vis":
        x_unit = 'Wavelength / nm'
        y_unit = '-log($R$) / a.u.'
    
    text = 'Path of your data'
    name_data = FileOpen(text)
    
    if '_PSD' in name_data and not '_iPW' in name_data:
        
        # data set (time values in first row; energy values in first column; each column is one spectrum)
        data_psd = pd.read_csv(r''+name_data, sep="\t", header=0)
        data_psd = data_psd.values
        
        phi = data_psd[0,1:] # time values [s]
        phi = np.reshape(phi,(phi.size,1))
        phi = phi.astype(int)
        
        # plot the PSD spectra
        plt.figure(figsize=(10,5))
        plt.plot(data_psd[1:,0],data_psd[1:,1:])
        plt.xlabel(x_unit) # get x axis label
        plt.ylabel(y_unit) # get y axis label
        #plt.ylim(np.amin(data_psd[1:,1:]), np.amax(data_psd[1:,1:]))
        plt.xlim(np.amin(data_psd[1:,0]), np.amax(data_psd[1:,0]))
        plt.legend(phi[:,0], title = r'$\varphi$ / °', loc = 'upper left', bbox_to_anchor=(1,1)) # legend outside frame
        
    elif '_averaged' in name_data and not 'course' in name_data:
        
        # Data of catalyst spectra
        data = pd.read_csv(r''+name_data, sep="\t", header = 0)
        data = data.values
        
        t_inp = data[0,1:] # time values [s]
        t_inp = np.reshape(t_inp,(t_inp.size,1))
        
        # plot the averaged spectra
        plt.figure(figsize=(10,5))
        plt.xlabel(x_unit) # get x axis label
        plt.ylabel(y_unit) # get y axis label
        #plt.ylim(np.amin(data[1:,1:]), np.amax(data[1:,1:]))
        plt.xlim(np.amin(data[1:,0]), np.amax(data[1:,0]))
        
        n_pl = 12 # maximum number of spectra to plot
        
        if n_pl < len(t_inp):
            t_legend = np.zeros((n_pl,1))
            i_pl = len(t_inp) // n_pl
            for i in range(n_pl):
                t_legend[i,0] = t_inp[(i+1)*i_pl-1,0]
                plt.plot(data[1:,0],data[1:,(i+1)*i_pl-1])
            plt.legend(t_legend[:n_pl,0], title = r'$t$ / s', loc = 'upper left', bbox_to_anchor=(1,1)) # legend outside frame
        
        else:
            plt.plot(data[1:,0],data[1:,1:])
            plt.legend(t_inp[:len(t_inp),0], title = r'$t$ / s', loc = 'upper left', bbox_to_anchor=(1,1)) # legend outside frame
        
    
    elif '_course' in name_data:
        
        n_sp = int(Entry_n_sp.get()) # number of spectra per period

        # Data of catalyst spectra
        data = pd.read_csv(r''+name_data, sep="\t", header = None)
        data = data.values

        t_inp = data[1:,0]
        t_inp = np.reshape(t_inp,(t_inp.size,1))
        t_inp = t_inp.astype(float)

        points = data[0,1:]
        points = np.reshape(points,(points.size,1))
        points = points.astype(float)

        # Plot the course
        for i in np.arange(0,len(points)):
            if i%12 == 0: # opens a new plot window every XXX lines (insert number of your choice). Otherwise the colours get confusing
                # Highlight the different phases of your periodic stimulation
                plt.figure(figsize=(10,5))
                for j in np.arange(min(t_inp),max(t_inp),n_sp*t_inp[0,0]): 
                    plt.axvspan(j, j+n_sp/2*t_inp[0,0], facecolor='k', alpha=0.25)
            
    elif '_iPW' in name_data:
            
            data = pd.read_csv(r''+name_data, sep="\t")#, header=0)
            data = data.values
            
            # plot the time constants
            plt.figure(figsize=(10,5))
            
            if len(data) > 100: # if more than xxx time constants plot as connected lines
                plt.plot(data[:,0], data[:,3])
            else: # else plot as dots
                plt.plot(data[:,0], data[:,3], 'o')
            plt.xlabel(x_unit) # get x axis label
            plt.ylabel('$t$ / s') # get y axis label
            
            
    return


def course():
    
    # x and y unit of the spectra
    if spectro.get() == "DRIFTS":
        x_unit = 'Wavenumber / 1/cm'
        y_unit = '-log($R$) / a.u.'
    elif spectro.get() == "Raman":
        x_unit = 'Raman Shift / 1/cm'
        y_unit = 'Raman Intensity / a.u.'
    elif spectro.get() == "UV-Vis":
        x_unit = 'Wavelength / nm'
        y_unit = '-log($R$) / a.u.'
    
    # input via the GUI 
    n_sp = int(Entry_n_sp.get()) # number of spectra per period
    n_per = int(Entry_n_per.get()) # number of periods that will be processed
    cutoff_per = int(Entry_cutoff_per.get()) # number of periods to cut off
    
    text = 'Path of your catalyst spectra' # choose data set
    name_data = FileOpen(text)
    
    text = 'Path of your points'
    name_points = FileOpen(text)
    
    points = np.genfromtxt(r''+name_points, delimiter="\n") # point positions
    points = np.sort(points) # sort points in ascending order
    
    
    if spectro.get() == "DRIFTS" and n_per != 1:
        
        # data set (wavenumber in first column; each column is one spectrum)
        data = pd.read_csv(r''+name_data, sep="\t", header = None)
        data = data.values
        
        # pick time values automatically 
        name_t = name_data.split('.')
        name_t = name_t[0] + '_t.' + name_t[1]
        
        t_inp = np.genfromtxt(r''+name_t, delimiter="\t") # time values [s]
        t_inp = np.delete(t_inp, np.s_[0,1,3,4], axis=1)
        
        energy_values = data[:,0] # wavenumbers [1/cm]    
        energy_values = np.reshape(energy_values,(energy_values.size,1)) # make 2D array for further computations
    
        data = np.delete(data, 0, axis=1) # delete energy values from data
        

    elif spectro.get() == "Raman" and n_per != 1:
        
        # data set (time values in first column; Raman shift in first row; each row is one spectrum)
        data = pd.read_csv(r''+name_data, sep="\t", header=None)
        data = data.values
        
        t_inp = data[1:,0] # time values [0-1]
        t_inp = np.reshape(t_inp,(t_inp.size,1))
        t_inp = (t_inp - t_inp[0]) * 24 * 3600 # time values [s]
        
        t_dummy = t_inp[1,0] + t_inp[len(t_inp)-1,0] # total experiment time [s]
        t_dummy = t_dummy/(n_sp*n_per) # averaged aquisition time [s]
        for i in range(0,len(t_inp)):
            t_inp[i,0] = (i+1)*t_dummy

        energy_values = data[0,1:] # Raman shifts [1/cm]
        energy_values = np.reshape(energy_values,(energy_values.size,1))
        
        data = np.delete(data, 0, axis=1) # delete time values from data
        data = np.delete(data, 0, axis=0) # delete energy values from data
        data = data.T
        
        # remove cosmic ray spikes
        if Entry_spike_rem.get() == 1:
            spike_factor = 1.2 # spikes with x times higher intensity than in spectrum before or after are removed
            for i in range(1,len(data.T)-1): # from second to second last spectrum 
                for j in range(0,len(data)):
                    if data[j,i] > spike_factor*data[j,i-1] and data[j,i] > spike_factor*data[j,i+1]:
                        data[j,i] = (data[j,i-1]+data[j,i+1])/2

        
    elif spectro.get() == "UV-Vis" and n_per != 1:
        
        # data set (time values in first row; wavelength in first column; each column is one spectrum)
        data = pd.read_csv(r''+name_data, sep=";", header=7)
        data = data.replace(',','.',regex=True)
        data = data.values

        t_inp = data[0,3:]  # time values [s]
        t_inp = np.reshape(t_inp,(t_inp.size,1))
        t_inp = t_inp.astype(float)
        t_inp[:,0] = t_inp[:,0] + t_inp[1,0]
        t_inp = t_inp/100000

        energy_values = data[1:,0] # wavelengths [nm]
        energy_values = np.reshape(energy_values,(energy_values.size,1))
        energy_values = energy_values.astype(float)

        data = np.delete(data, [0,1,2], axis=1) # delete energy values from data
        data = np.delete(data, 0, axis=0) # delete time values from data
        data = data.astype(float)
        
    elif n_per == 1:
        
        # data set (time values in first row; energy values in first row; each column is one spectrum)
        data = pd.read_csv(r''+name_data, sep="\t", header=0)
        data = data.values

        t_inp = data[0,1:] # time values [s]
        t_inp = np.reshape(t_inp,(t_inp.size,1))

        energy_values = data[1:,0] # Raman shift [1/cm]
        energy_values = np.reshape(energy_values,(energy_values.size,1))

        data = np.delete(data, 0, axis=1) # delete energy values from data
        data = np.delete(data, 0, axis=0) # delete time values from data
    
    
    # only for data that were not averaged before
    if n_per != 1:
        # cut off spectra at the end
        if n_per != len(data.T)/n_sp:
            cutoff_sp = (int(len(data.T)/n_sp) - n_per) * n_sp
            data = np.delete(data, np.s_[-(cutoff_sp):], axis=1)
            t_inp = np.delete(t_inp,np.s_[-(cutoff_sp):], axis=0)
        
        # cut off spectra from the beginning
        if cutoff_per != 0:
            cutoff_sp = cutoff_per * n_sp # Calculated number of spectra to cut off
            data = np.delete(data, np.s_[:cutoff_sp], axis=1)
            t_inp = np.delete(t_inp,np.s_[-(cutoff_sp):], axis=0)
          
        # bring energy_values in an ascending order and sort spectra correctly
        if energy_values[0,0] > energy_values[1,0]:
            energy_values = np.flip(energy_values)
            data = np.flip(data, axis=0)

    # normalize all averaged spectra to the single highest peak
    if spectro.get() == "Raman" and Entry_norm.get() == 1:
        I_max = np.max(data)
        # print("I_max = " + str(I_max)) # print the highest intensity to normalize other ranges
        data = data/I_max
            
    
    i = 0
    for val in points:
        points[i] = min(energy_values, key=lambda x:abs(x-val)) #If the wavenumber from 'points' and data[:,0] don't fit 100% the value with the lowest deviation will be taken
        i = i+1
        
    # Plot the course
    pos = np.zeros(len(points))
    for i in np.arange(0,len(points)):
        if i%12 == 0: # opens a new plot window every XXX lines (insert number of your choice). Otherwise the colours get confusing
            # highlights the different phases of your periodic stimulation
            plt.figure(figsize=(10,5))
            for j in np.arange(min(t_inp/60),max(t_inp/60),n_sp*t_inp[0,0]/60): #divide by 60 to get from s to min
                plt.axvspan(j, j+n_sp/2*t_inp[0,0]/60, facecolor='k', alpha=0.25)
        
        dummy = np.where(energy_values == points[i])
        pos[i] = dummy[0]
        
        plt.plot(t_inp/60,data[int(pos[i]),:], label = str(int(np.around(points[i],0))))
        plt.legend(title = x_unit, loc='upper right')
        plt.xlabel('$t$ / min')
        plt.ylabel(y_unit) # Gets y axis label
        plt.xlim(np.amin(t_inp/60), np.amax(t_inp/60))

    # Write sth. to save output as txt
    output1 = pd.DataFrame({'t / min': t_inp[:,0]/60})

    output2 = pd.DataFrame(data=data[pos.astype(int),:].T, columns=np.around(points,0).astype(int))

    output = pd.merge(output1,output2, left_index=True, right_index=True)

    # Save data frame to a file?

    text = 'Shall the temporal courses be saved as .txt?'
    name = name_data.split('.') #Will be used as filename
    name = name[0] + '_course.txt'
    yesno(name, output, text)
    
    return

def time_resolved(): # averaging of TRS
    
    # input via the GUI 
    n_sp = int(Entry_n_sp.get()) # number of spectra per period
    n_per = int(Entry_n_per.get()) # number of periods that will be processed
    cutoff_per = int(Entry_cutoff_per.get()) # number of periods to cut off
    
    text = 'Path of your catalyst spectra' # choose data set
    name_data = FileOpen(text)
    
    start = time.time()
    
    
    if spectro.get() == "DRIFTS":
        
        # data set (wavenumber in first column; each column is one spectrum)
        data = pd.read_csv(r''+name_data, sep="\t", header = None)
        data = data.values
        
        # pick time values automatically
        name_t = name_data.split('.')
        name_t = name_t[0] + '_t.' + name_t[1]
        
        t_inp = np.genfromtxt(r''+name_t, delimiter="\t") # time values [s]
        t_inp = np.delete(t_inp, np.s_[0,1,3,4], axis=1)
        
        energy_values = data[:,0] # wavenumbers [1/cm]    
        energy_values = np.reshape(energy_values,(energy_values.size,1)) # make 2D array for further computations
    
        data = np.delete(data, 0, axis=1) # delete energy values from data
        

    elif spectro.get() == "Raman":
        
        # data set (time values in first column; Raman shift in first row; each row is one spectrum)
        data = pd.read_csv(r''+name_data, sep="\t", header=None)
        data = data.values
        
        t_inp = data[1:,0] # time values [0-1]
        t_inp = np.reshape(t_inp,(t_inp.size,1))
        t_inp = (t_inp - t_inp[0]) * 24 * 3600 # time values [s]
        
        t_dummy = t_inp[1,0] + t_inp[len(t_inp)-1,0] # total experiment time [s]
        t_dummy = t_dummy/(n_sp*n_per) # averaged aquisition time [s]
        for i in range(0,len(t_inp)):
            t_inp[i,0] = (i+1)*t_dummy

        energy_values = data[0,1:] # Raman shifts [1/cm]
        energy_values = np.reshape(energy_values,(energy_values.size,1))
        
        data = np.delete(data, 0, axis=1) # delete time values from data
        data = np.delete(data, 0, axis=0) # delete energy values from data
        data = data.T
        
        # remove cosmic ray spikes
        if Entry_spike_rem.get() == 1:
            spike_factor = 1.2 # spikes with x times higher intensity than in spectrum before or after are removed
            for i in range(1,len(data.T)-1): # from second to second last spectrum 
                for j in range(0,len(data)):
                    if data[j,i] > spike_factor*data[j,i-1] and data[j,i] > spike_factor*data[j,i+1]:
                        data[j,i] = (data[j,i-1]+data[j,i+1])/2

        
    elif spectro.get() == "UV-Vis":
        
        # data set (time values in first row; wavelength in first column; each column is one spectrum)
        data = pd.read_csv(r''+name_data, sep=";", header=7)
        data = data.replace(',','.',regex=True)
        data = data.values

        t_inp = data[0,3:]  # time values [s]
        t_inp = np.reshape(t_inp,(t_inp.size,1))
        t_inp = t_inp.astype(float)
        t_inp[:,0] = t_inp[:,0] + t_inp[1,0]
        t_inp = t_inp/100000

        energy_values = data[1:,0] # wavelengths [nm]
        energy_values = np.reshape(energy_values,(energy_values.size,1))
        energy_values = energy_values.astype(float)

        data = np.delete(data, [0,1,2], axis=1) # delete energy values from data
        data = np.delete(data, 0, axis=0) # delete time values from data
        data = data.astype(float)
        
   
    # cut off spectra at the end
    if n_per != len(data.T)/n_sp:
        cutoff_sp = (int(len(data.T)/n_sp) - n_per) * n_sp
        data = np.delete(data, np.s_[-(cutoff_sp):], axis=1)
        t_inp = np.delete(t_inp,np.s_[-(cutoff_sp):], axis=0)
    
    # cut off spectra from the beginning
    if cutoff_per != 0:
        cutoff_sp = cutoff_per * n_sp # Calculated number of spectra to cut off
        data = np.delete(data, np.s_[:cutoff_sp], axis=1)
        t_inp = np.delete(t_inp,np.s_[-(cutoff_sp):], axis=0)
          
    # average all periods into one period
    data = np.split(data, n_per-cutoff_per, axis=1) # split the wholeness of all spectra into minor ndarrays for each period    
    data = np.divide(sum(data),n_per) # sum up all cells of the created ndarrays that have the same index and divide by the number of periods       
        
    # bring energy_values in an ascending order and sort spectra correctly
    if energy_values[0,0] > energy_values[1,0]:
        energy_values = np.flip(energy_values)
        data = np.flip(data, axis=0)

    # normalize all averaged spectra to the single highest peak
    if spectro.get() == "Raman" and Entry_norm.get() == 1:
        I_max = np.max(data)
        # print("I_max = " + str(I_max)) # print the highest intensity to normalize other ranges
        data = data/I_max
    
    data = np.concatenate((energy_values, data), axis=1) # concatenate energy values and averaged spectra

    t_avg = np.insert(t_inp[:n_sp,0], 0, 0) # insert 0 for later saving
    t_avg = np.reshape(t_avg,(t_avg.size,1))
    data = np.concatenate((t_avg.T, data), axis=0) # concatenate averaged spectra and time values
    
    end = time.time()
    print('Runtime: ' + str(end-start) +' s.')
    
    
    # x and y unit of the spectra
    if spectro.get() == "DRIFTS":
        x_unit = 'Wavenumber / 1/cm'
        y_unit = '-log($R$) / a.u.'
    elif spectro.get() == "Raman":
        x_unit = 'Raman Shift / 1/cm'
        y_unit = 'Raman Intensity / a.u.'
    elif spectro.get() == "UV-Vis":
        x_unit = 'Wavelength / nm'
        y_unit = '-log($R$) / a.u.'
    
    
    # save the averaged spectra
    t_inp = np.around(t_inp,1)
    t_str = [' (' + str(item) + ' s)' for item in t_inp[:n_sp,0]]
    for i in range(0,n_sp):
        t_str[i] = '#' + str(i+1) + t_str[i]
    header = [x_unit] + t_str
    output = pd.DataFrame(data, columns = header)

    text = 'Shall the averaged spectra be saved as .txt?'
    name = name_data.split('.')
    name = name[0] + '_' + str(cutoff_per) + '_periods_cutoff_' + 'averaged.txt'
    yesno(name, output, text)

    # plot the averaged spectra
    plt.figure(figsize=(10,5))
    plt.xlabel(x_unit) # get x axis label
    plt.ylabel(y_unit) # get y axis label
    #plt.ylim(np.amin(data[1:,1:]), np.amax(data[1:,1:]))
    plt.xlim(np.amin(data[1:,0]), np.amax(data[1:,0]))
    
    n_pl = 12 # maximum number of spectra to plot
    
    if n_pl < n_sp:
        t_legend = np.zeros((n_pl,1))
        i_pl = n_sp // n_pl
        for i in range(n_pl):
            t_legend[i,0] = t_inp[(i+1)*i_pl-1,0]
            plt.plot(data[1:,0],data[1:,(i+1)*i_pl-1])
        plt.legend(t_legend[:n_pl,0], title = r'$t$ / s', loc = 'upper left', bbox_to_anchor=(1,1)) # legend outside frame
    
    else:
        plt.plot(data[1:,0],data[1:,1:])
        plt.legend(t_inp[:n_sp,0], title = r'$t$ / s', loc = 'upper left', bbox_to_anchor=(1,1)) # legend outside frame
    
    
    return

def Baseline():
    
    text = 'Path of your averaged spectra'
    name_data = FileOpen(text)
    
    text = 'Path of your baseline points'
    name_points = FileOpen(text)
    
    start = time.time()

    # Data of catalyst spectra
    data = pd.read_csv(r''+name_data, sep="\t", header = 0)
    data = data.values

    t_inp = data[0,1:] # time values [s]
    t_inp = np.reshape(t_inp,(t_inp.size,1))

    energy_values = data[1:,0]
    energy_values = np.reshape(energy_values,(energy_values.size,1))

    data = np.delete(data, 0, axis=1) # delete energy values from data
    data = np.delete(data, 0, axis=0) # delete time values from data

    # generating baseline
    point_pos = np.genfromtxt(r''+name_points, delimiter="\n")

    # Compares every value in peak_pos with the wavenumbers from psd_spectra and the closest value is taken
    index_E = np.zeros(len(point_pos))
    i = 0
    for val in point_pos:      
        index_E[i] = abs(energy_values-point_pos[i]).argmin()
        i = i+1

    I_dummy = np.zeros((len(t_inp), len(point_pos)+2))
    I_bl = np.zeros((len(energy_values), len(t_inp)))

    for i in range(len(point_pos)+1):
        if i == 0:
            I_dummy[:,i] = data[0,:]
            I_dummy[:,i+1] = data[int(index_E[i]),:]
            m = (I_dummy[:,i] - I_dummy[:,i+1]) / (energy_values[0] - energy_values[int(index_E[i])])
            b = I_dummy[:,i] - m * energy_values[0]
            I_bl[:int(index_E[i]),:] = m * energy_values[:int(index_E[i])] + b
        elif i > 0 and i < len(point_pos):
            I_dummy[:,i+1] = data[int(index_E[i]),:]
            m = (I_dummy[:,i+1] - I_dummy[:,i]) / (energy_values[int(index_E[i])] - energy_values[int(index_E[i-1])])
            b = I_dummy[:,i+1] - m * energy_values[int(index_E[i])]
            I_bl[int(index_E[i-1]):int(index_E[i]),:] = m * energy_values[int(index_E[i-1]):int(index_E[i])] + b
        else:
            I_dummy[:,i+1] = data[len(energy_values)-1,:]
            m = (I_dummy[:,i+1] - I_dummy[:,i]) / (energy_values[len(energy_values)-1] - energy_values[int(index_E[i-1])])
            b = I_dummy[:,i] - m * energy_values[int(index_E[i-1])]
            I_bl[int(index_E[i-1]):,:] = m * energy_values[int(index_E[i-1]):] + b
    
    
    I_bl = np.concatenate((energy_values, I_bl), axis=1) # concatenate energy values and averaged spectra

    t_avg = np.insert(t_inp, 0, 0) # insert 0 for later saving
    t_avg = np.reshape(t_avg,(t_avg.size,1))
    I_bl = np.concatenate((t_avg.T, I_bl), axis=0) # concatenate averaged spectra and time values
    
    end = time.time()
    print('Runtime: ' + str(end-start) +' s.')
    
    
    # x and y unit of the spectra
    if spectro.get() == "DRIFTS":
        x_unit = 'Wavenumber / 1/cm'
        y_unit = '-log($R$) / a.u.'
    elif spectro.get() == "Raman":
        x_unit = 'Raman Shift / 1/cm'
        y_unit = 'Raman Intensity / a.u.'
    elif spectro.get() == "UV-Vis":
        x_unit = 'Wavelength / nm'
        y_unit = '-log($R$) / a.u.'
        
    # save the baselines
    t_str = np.around(t_inp[:,0], 1)
    t_str = [' (' + str(item) + ' s)' for item in t_str]
    for i in range(0,t_inp.size):
        t_str[i] = '#' + str(i+1) + t_str[i]
    header = [x_unit] + t_str
    output = pd.DataFrame(I_bl, columns = header)

    text = 'Shall the baselines be saved as .txt?'
    name = name_data.split('.')
    name = name[0] + '_baseline.txt'
    yesno(name, output, text)
       
    # plot the the last spectrum of the first gasphase and the corresponding baseline and difference spectrum
    plt.figure(figsize=(10,5))
    plt.xlabel(x_unit) # get x axis label
    plt.ylabel(y_unit) # get y axis label
    #plt.ylim(np.amin(data[1:,1:]), np.amax(data[1:,1:]))
    plt.xlim(np.amin(energy_values), np.amax(energy_values))

    plt.plot(energy_values, data[:,int(len(t_inp)/2)-1])
    plt.plot(energy_values, I_bl[1:,int(len(t_inp)/2)])
    plt.plot(energy_values, data[:,int(len(t_inp)/2)-1] - I_bl[1:,int(len(t_inp)/2)])

    # plt.legend(['spectrum', 'baseline', 'difference'], loc = 'upper left', bbox_to_anchor=(1,1)) # legend outside frame
    plt.legend(['spectrum', 'baseline', 'difference'], loc = 'upper right') # legend inside frame
    
    
    return


'''
_______________________________________________________________________________
GUI stuff
_______________________________________________________________________________
'''

PSD_GUI = Tk()
PSD_GUI.title('Have fun with PSD!')

frame_left = Frame(PSD_GUI)
frame_left.pack(side=LEFT)

frame_right = Frame(PSD_GUI)
frame_right.pack(side=RIGHT)




Label_DD_spectro = Label(frame_left, text = 'Choose the spectroscopy type!').pack()
spectro_list = ['DRIFTS', 'Raman', 'UV-Vis']
spectro = StringVar(PSD_GUI)
spectro.set('DRIFTS')
DD_spectro = OptionMenu(frame_left, spectro, *spectro_list)
DD_spectro.pack()

Label_n_sp = Label(frame_left, text = 'Type in the number of spectra per period!').pack()
Entry_n_sp = StringVar()
Entry_n_sp = Entry(frame_left, textvariable = Entry_n_sp)
Entry_n_sp.insert(END,'0')
Entry_n_sp.pack()

Label_n_per = Label(frame_left, text = 'Type in the number of periods!').pack()
Entry_n_per = StringVar()
Entry_n_per = Entry(frame_left, textvariable = Entry_n_per)
Entry_n_per.insert(END,'0')
Entry_n_per.pack()

Label_cutoff_per = Label(frame_left, text = 'Choose the number of periods to cut off!').pack()
Entry_cutoff_per = StringVar()
Entry_cutoff_per = Entry(frame_left, textvariable = Entry_cutoff_per)
Entry_cutoff_per.insert(END,'0')
Entry_cutoff_per.pack()

Label_dphi = Label(frame_left, text = 'Choose your phase resolution!').pack()
Entry_dphi = StringVar()
Entry_dphi = Entry(frame_left, textvariable = Entry_dphi)
Entry_dphi.insert(END,'30')
Entry_dphi.pack()

Label_k_harmonic = Label(frame_left, text = 'Type in the harmonic to demodulate with \n (1, 3, 5, ... for sine or 0 for rectangular function)!').pack()
Entry_k_harmonic = StringVar()
Entry_k_harmonic = Entry(frame_left, textvariable = Entry_k_harmonic)
Entry_k_harmonic.insert(END,'1')
Entry_k_harmonic.pack()

Label_SpikeRemoving = Label(frame_left, text = 'Only relevant for Raman spectroscopy: \n Do you want to remove cosmic ray spikes?').pack()
Entry_spike_rem = IntVar(PSD_GUI)
Entry_spike_rem.set(value=1)
Box_SpikeRemoving = Checkbutton(frame_left, text = 'spike removing', variable=Entry_spike_rem, onvalue=1, offvalue=0)
Box_SpikeRemoving.pack()

Label_Normalization = Label(frame_left, text = 'Do you want to normalize the data?').pack()
Entry_norm = IntVar(PSD_GUI)
Entry_norm.set(value=1)
Box_Normalization = Checkbutton(frame_left, text = 'normalization', variable=Entry_norm, onvalue=1, offvalue=0)
Box_Normalization.pack()

Label_PSD_calc = Label(frame_right, text = 'Calculate PSD spectra from time-resolved spectra (TRS):').pack()
Bt_PSD_calc = Button(frame_right, text = 'PSD', command = PSD_calc).pack()

Label_time_resolved = Label(frame_right, text = 'Average the TRS into one period without PSD:').pack()
Bt_time_resolved = Button(frame_right, text = 'TRS Averaging', command = time_resolved).pack()

Label_PointPicking = Label(frame_right, text = 'If a data point is clicked on, all data points \n clicked on until now are written into a file \n (If no graph is shown, resize the window!):').pack()
Bt_PointPicking = Button(frame_right, text = 'point picking', command = PointPicking).pack()

Label_in_phase_angle = Label(frame_right, text = 'Calculate in-phase angle and in-phase time:').pack()
Bt_in_phase_angle = Button(frame_right, text = 'in phase angle', command = in_phase_angle).pack()

Label_Spectra_diff = Label(frame_right, text = 'Here you calculate difference spectra:').pack()
Bt_Spectra_diff = Button(frame_right, text = 'difference spectra', command = Spectra_diff).pack()

Label_Baseline = Label(frame_right, text = 'Here you gernerate baselines:').pack()
Bt_Baseline = Button(frame_right, text = 'baseline', command = Baseline).pack()

Label_Show_Graph = Label(frame_right, text = 'Plot graphs you like:').pack()
Bt_Show_Graph = Button(frame_right, text = 'show graph', command = Show_Graph).pack()

Label_Show_Points = Label(frame_right, text = 'Here you can highlight chosen point positions:').pack()
Bt_Show_Points = Button(frame_right, text = 'show points', command = Show_Points).pack()

Label_course = Label(frame_right, text = 'Create course plots of chosen spectra at chosen point positions:').pack()
Bt_course = Button(frame_right, text = 'course plot', command = course).pack()






# Label_DD_spectro = Label(PSD_GUI, text = 'Choose the spectroscopy type!').pack()
# spectro_list = ['DRIFTS', 'Raman', 'UV-Vis']
# spectro = StringVar(PSD_GUI)
# spectro.set('DRIFTS')
# DD_spectro = OptionMenu(PSD_GUI, spectro, *spectro_list)
# DD_spectro.pack()

# Label_n_sp = Label(PSD_GUI, text = 'Type in the number of spectra per period!').pack()
# Entry_n_sp = StringVar()
# Entry_n_sp = Entry(PSD_GUI, textvariable = Entry_n_sp)
# Entry_n_sp.insert(END,'0')
# Entry_n_sp.pack()

# Label_n_per = Label(PSD_GUI, text = 'Type in the number of periods!').pack()
# Entry_n_per = StringVar()
# Entry_n_per = Entry(PSD_GUI, textvariable = Entry_n_per)
# Entry_n_per.insert(END,'0')
# Entry_n_per.pack()

# Label_cutoff_per = Label(PSD_GUI, text = 'Choose the number of periods to cut off!').pack()
# Entry_cutoff_per = StringVar()
# Entry_cutoff_per = Entry(PSD_GUI, textvariable = Entry_cutoff_per)
# Entry_cutoff_per.insert(END,'0')
# Entry_cutoff_per.pack()

# Label_dphi = Label(PSD_GUI, text = 'Choose your phase resolution!').pack()
# Entry_dphi = StringVar()
# Entry_dphi = Entry(PSD_GUI, textvariable = Entry_dphi)
# Entry_dphi.insert(END,'30')
# Entry_dphi.pack()

# Label_k_harmonic = Label(PSD_GUI, text = 'Type in the harmonic to demodulate with \n (1, 3, 5, ... for sine or 0 for rectangular function)!').pack()
# Entry_k_harmonic = StringVar()
# Entry_k_harmonic = Entry(PSD_GUI, textvariable = Entry_k_harmonic)
# Entry_k_harmonic.insert(END,'1')
# Entry_k_harmonic.pack()

# Label_SpikeRemoving = Label(PSD_GUI, text = 'Do you want to remove cosmic ray spikes?').pack()
# Entry_spike_rem = IntVar(PSD_GUI)
# Entry_spike_rem.set(value=1)
# Box_SpikeRemoving = Checkbutton(PSD_GUI, text = 'spike removing', variable=Entry_spike_rem, onvalue=1, offvalue=0)
# Box_SpikeRemoving.pack()

# Label_Normalization = Label(PSD_GUI, text = 'Do you want to normalize the data?').pack()
# Entry_norm = IntVar(PSD_GUI)
# Entry_norm.set(value=1)
# Box_Normalization = Checkbutton(PSD_GUI, text = 'normalization', variable=Entry_norm, onvalue=1, offvalue=0)
# Box_Normalization.pack()

# Label_PSD_calc = Label(PSD_GUI, text = 'Calculate PSD spectra from time resolved ones:').pack()
# Bt_PSD_calc = Button(PSD_GUI, text = 'PSD', command = PSD_calc).pack()

# Label_time_resolved = Label(PSD_GUI, text = 'Average the TRS into one period:').pack()
# Bt_time_resolved = Button(PSD_GUI, text = 'TRS Averaging', command = time_resolved).pack()

# Label_PointPicking = Label(PSD_GUI, text = 'If no graph is shown, resize the window! \n If a data point is clicked on, all data points \n clicked on until now are written into a file:').pack()
# Bt_PointPicking = Button(PSD_GUI, text = 'point picking', command = PointPicking).pack()

# Label_in_phase_angle = Label(PSD_GUI, text = 'Calculate in phase angle and in phase time:').pack()
# Bt_in_phase_angle = Button(PSD_GUI, text = 'in phase angle', command = in_phase_angle).pack()

# Label_Spectra_diff = Label(PSD_GUI, text = 'Here you calculate difference spectra:').pack()
# Bt_Spectra_diff = Button(PSD_GUI, text = 'difference spectra', command = Spectra_diff).pack()

# Label_Baseline = Label(PSD_GUI, text = 'Here you gernerate baselines:').pack()
# Bt_Baseline = Button(PSD_GUI, text = 'baseline', command = Baseline).pack()

# Label_Show_Graph = Label(PSD_GUI, text = 'Plot graphs you like:').pack()
# Bt_Show_Graph = Button(PSD_GUI, text = 'show graph', command = Show_Graph).pack()

# Label_Show_Points = Label(PSD_GUI, text = 'here you can highlight chosen point positions:').pack()
# Bt_Show_Points = Button(PSD_GUI, text = 'show points', command = Show_Points).pack()

# Label_course = Label(PSD_GUI, text = 'Create course plots of chosen spectra at chosen point positions:').pack()
# Bt_course = Button(PSD_GUI, text = 'course plot', command = course).pack()






PSD_GUI.mainloop()