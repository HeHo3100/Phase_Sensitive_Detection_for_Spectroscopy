# Phase Sensitive Detection for Spectroscopy

https://doi.org/10.5281/zenodo.3613876

A GUI with multiple features that executes Phase Sensitive Detection in Modulation Excitation Spectroscopy to get from time
resolved data to phase resolved spectra.
The raw data used by me were collected via DRIFT Spectroscopy but this script is applicable to all other time resolved
spectroscopic techniques.

A spectroscopic data processing tool by Jakob Weyel, Eduard-Zintl-Institut für Anorganische und Physikalische Chemie,
TU Darmstadt

For further questions feel free to contact me via mail: jakob.weyel@tu-darmstadt.de

For 'plug and play' use install anaconda as your python and use scipy as I do.

Executing this script opens a GUI with several buttons:

#### TRS->PSD:
Reads time resolved spectral data, processes them according to PSD and saves the resulting phase resolved spectra as .txt.
Every column of the raw data has to contain the intensity of one spectrum, the first column contains e. g. the respective
energy/frequency/wavenumber/wavelength to do the processing correctly. Needs the correct number of spectra collected during
one modulation period to process your data correctly (first textbox of the GUI)! The second textbox lets you decide how many
modulation periods from the beginning are cut of. The third textbox lets you type in the desired phase resolution of your
Output spectra in °.
###### Input needed:
RawData.txt, ReferenceData.txt of the same size (only if you want to process modulated difference spectra, otherwise press Esc), time values as a file with the same name as the raw data with the addition of _t ('FilenameOfYourRawData_t.txt') which contains temporal data as acolumn vector, correct 'number of spectra per period' in first textbox
###### Output:
FilenameOfYourRawData_PSD.txt
  
#### difference spectra:
Import two data sets of the same size which contain spectra of your choice and calculate their difference (Spectrum1-Spectrum2)
###### Input needed:
RawData.txt, ReferenceData.txt of the same size
###### Output:
FilenameOfYourRawData_TRS_Spektrum.txt
  
#### peak picking:
Imports a spectrum of your choice in which you can click on points of your choice in the spectrum to save all chosen
energy/frequency/wavenumber/wavelength into a.txt file. After clicking for the first time and thus creating your file,
every other click overwrites the existing file with the data points which were clicked on in the current window.
###### Input needed:
RawData.txt or ProcessedData.txt (any kind of spectrum suffices)
###### Output:
FilenameOfYourRawData_Banden.txt
  
#### show peaks:
Import a data set with energy/frequency/wavenumber/wavelength values (which you may create e. g. via the button 'peak
picking') and highlight these positions in your latest matplotlib window as red vertical lines.
###### Input needed:
AnyData.txt, AnyData_Banden.txt
###### Output:
None
  
#### in phase angle:
Import a data set of PSD spectra, peak positions of your choice and the temporal data of your raw data related to your PSD
spectra to determine for every band at which phase angle it has its maximum value. Temporal data and a correct input in the
first textbox (number of spectra per period) is needed to calculate time values out of the maximum phase angles.
###### Input needed:
PSDData.txt, AnyData_Banden.txt, RawData_t.txt, correct 'number of spectra per period' in first textbox
###### Output:
PSDData_Banden_iPW.txt
  
#### show graph:
Plots a data set of your choice.
###### Input needed:
AnyData.txt
###### Output:
None

#### contour plot:
Plots the temporal course of one whole phase resolved spectra set and links the appearance of bands with the time of one period.
Thus the correct time vector of your raw data is needed as an ascii file as well as the correct 'number of spectra per period'
as input for the first textbox
###### Input needed:
PSDData.txt, time values of your raw data (e. g. 'FilenameOfYourRawData_t.txt') which contains temporal data as acolumn vector, correct 'number of spectra per period' in first textbox
###### Output:
None
  
#### course plot:
Plots the temporal course of bands of your choice from your time resolved data in 2D I vs t diagram.
###### Input needed:
RawData.txt, ReferenceData.txt of the same size (only if you want to process modulated difference spectra, otherwise press Esc), time values as a file with the same name as the raw data with the addition of _t ('FilenameOfYourRawData_t.txt') which contains temporal data as acolumn vector and AnyData_Banden.txt, 
###### Output:
None

