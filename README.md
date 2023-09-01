https://doi.org/10.5281/zenodo.3613876

A GUI with multiple features that executes Phase Sensitive Detection in Modulation Excitation Spectroscopy to get from time
resolved data to phase resolved spectra.

A spectroscopic data processing tool by Jakob Weyel, Eduard-Zintl-Institut für Anorganische und Physikalische Chemie,
TU Darmstadt

The raw data used by me were collected doing DRIFT Spectroscopy but this script is applicable to all other time resolved
spectroscopic techniques, just make sure the format of your data fits the script: time resolved spectra should be stored in
one textfile with the first column being the energy/wavenumber/wavelength/etc. values and the spectra following column-wise.
Time values of your experiment have to be stored in a separate textfile containing a column vector.

When you load data into this script, make sure that the columns are separated by tab stops.

For 'plug and play' use install anaconda as your python and use spyder as I do.

For further questions feel free to contact me via mail: jakob.weyel@tu-darmstadt.de

Executing this script opens a GUI with several buttons:

#### TRS->PSD:
Reads time resolved spectral data, processes them according to PSD and saves the resulting phase resolved spectra as .txt.
Every column of the raw data has to contain the intensity of one spectrum, the first column contains e. g. the respective
energy/frequency/wavenumber/wavelength to do the processing correctly. Needs the harmonic you want to demodulate the spectra
with, number of periods, correct number of spectra collected during one modulation period, how many modulation periods from
the beginning are cut off and the desired phase resolution of your output spectra in ° to process your data correctly
(see textboxes of the GUI!!).
###### Input needed:
RawData.txt, ReferenceData.txt of the same size (only if you want to process modulated difference spectra, otherwise press Esc),
time values as a file with the same name as the raw data with the addition of _t ('FilenameOfYourRawData_t.txt') which contains
temporal data as a column vector, in textboxes: 'harmonic to demodulate with', 'number of periods', 'number of spectra per period',
'phase resolution'
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
FilenameOfYourRawData_peaks.txt
  
#### show peaks:
Import a data set with energy/frequency/wavenumber/wavelength values (which you may create e. g. via the button 'peak
picking') and highlight these positions in your latest matplotlib window as red vertical lines.
###### Input needed:
AnyData.txt, AnyData_peaks.txt
###### Output:
None
  
#### in phase angle:
Import a data set of PSD spectra, peak positions of your choice and the temporal data of your raw data related to your PSD
spectra to determine for every band at which phase angle it has its maximum value. Temporal data and a correct input in the
first textbox (number of spectra per period) is needed to calculate time values out of the maximum phase angles.
###### Input needed:
PSDData.txt, AnyData_peaks.txt, RawData_t.txt, in textboxes: 'number of spectra per period'
###### Output:
PSDData_peaks_iPW.txt
  
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
PSDData.txt, time values of your raw data (e. g. 'FilenameOfYourRawData_t.txt') which contains temporal data as a column vector,
correct 'number of spectra per period' in first textbox
###### Output:
None
  
#### course plot:
Plots the temporal course of bands of your choice from your time resolved data in 2D I vs t diagram.
###### Input needed:
RawData.txt, ReferenceData.txt of the same size (only if you want to process modulated difference spectra, otherwise hit Esc)
time values as a file with the same name as the raw data with the addition of _t ('FilenameOfYourRawData_t.txt') which contains
temporal data as acolumn vector and AnyData_peaks.txt, 
###### Output:
None
  
#### time resolved:
Creates a data set of time resolved spectra in a way that's easier to plot and to digest. All periods are averaged into one
'overall' period and only ten spectra of this period are taken.
###### Input needed:
RawData.txt, in textboxes: 'number of periods', 'number of spectra per period'
###### Output:
FilenameOfYourRawData_1period.txt
