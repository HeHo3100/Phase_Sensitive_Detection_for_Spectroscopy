https://doi.org/10.5281/zenodo.3613876

A GUI with multiple features that executes Phase Sensitive Detection (PSD) in Modulation Excitation Spectroscopy (MES) to get from
time-resolved data to phase-resolved spectra.

A spectroscopic data processing tool by Henrik Hoyer, Eduard-Zintl-Institut für Anorganische und Physikalische Chemie, TU Darmstadt
Originally by Jakob Weyel, formerly Eduard-Zintl-Institut für Anorganische und Physikalische Chemie, TU Darmstadt

The raw data used by me were collected doing DRIFT Spectroscopy (Bruker, OPUS software),
Raman spectroscopy (Horiba Jobin Yvon, LabSpec6 software) and UV-Vis spectroscopy (Avantes, AvaSoft software),
but this script is applicable to all other time-resolved spectroscopic techniques.
Just make sure the format of your data fits the script.

For 'plug and play' use install anaconda as your python and use spyder as I do. I recommend to set your
graphics backend to 'Automatic 'for functions such as 'PointPicking' and 'Show_Points' to work properly.

For further questions feel free to contact me via mail: henrik.hoyer@tu-darmstadt.de

Executing this script opens a GUI with several buttons:

#### PSD:
Reads time-resolved spectral data, processes them according to PSD and saves the resulting phase-resolved spectra as .txt.
Needs the correct spectroscopy type, number of spectra collected during one modulation period, number of periods you want to process
(choose 1 if the data are already averaged), how many modulation periods from the beginning are cut off, the desired phase resolution
of your output spectra in ° and the harmonic you want to demodulate the spectra with to process your data correctly
(see textboxes of the GUI!!). For Raman spectroscopy, choose if you want to have cosmic ray spikes removed and your data normalized.
The time-resolved raw data has to have the following format:
####### DRIFT spectroscopy: according to OPUS by Bruker
Every column of the raw data has to contain the intensity of one spectrum, the first column contains the respective wavenumber.
Tap stops as seperator. An extra file with time values in a column vector is required.
####### Raman spectroscopy: according to LabSpec6 by Horiba Jobin Yvon
Every column of the raw data has to contain the intensity of one spectrum, the first line contains the respective Raman shift
and the first column contains the time values. Tap stops as seperator.
####### UV-Vis spectroscopy: according to AvaSoft by Avantes
Every column of the raw data has to contain the intensity of one spectrum, the first column contains the respective wavelength
and the first line contains the time values. ';' as seperator.
####### Already averaged spectra (independebnt from sprectroscopy type):
Every column of the raw data has to contain the intensity of one spectrum, the first column contains the respective
wavenumber/Raman shift/wavelength and the first line contains the time values. Tap stops as seperator.
###### Input needed:
RawData.txt, (for DRIFTS: time values as a file with the same name as the raw data with the addition of _t ('FilenameOfYourRawData_t.txt')
which contains temporal data as a column vector), in textboxes: 'spectroscopy type', 'number of spectra per period', 'number of periods',
'number of periods to cut off', 'phase resolution', 'harmonic to demodulate with', for Raman: 'spike removing' and 'normalization'
###### Output:
FilenameOfYourRawData_PSD.txt

#### TRS Averaging:
Creates a data set of time-resolved spectra. All periods are averaged into one 'overall' period without doing PSD (the averaging is done
automatically in the 'PSD' function.
###### Input needed:
RawData.txt, in textboxes: 'spectroscopy type', 'number of spectra per period', 'number of periods', 'number of periods to cut off',
for Raman: 'spike removing' and 'normalization'
###### Output:
FilenameOfYourRawData_averaged.txt

#### point picking:
Imports a spectrum of your choice in which you can click on points of your choice in the spectrum to save all chosen
wavenumber/Raman shift/wavelength into a .txt file. After clicking for the first time and thus creating your file,
every other click overwrites the existing file with the data points which were clicked on in the current window.
###### Input needed:
RawData.txt or ProcessedData.txt (any kind of spectrum suffices), in textboxes: 'spectroscopy type'
###### Output:
FilenameOfYourRawData_points.txt

#### in-phase angle:
Import a data set of PSD spectra, peak positions of your choice and the averaged data related to your PSD for the time values to
determine for every band at which phase angle it has its maximum value, i.e., the in-hase angle.
###### Input needed:
PSDData.txt, AnyData_points.txt, AveragedData.txt, in textboxes: 'spectroscopy type'
###### Output:
PSDData_points_iPW.txt

#### difference spectra:
Import two averaged data sets of the same size which contain spectra of your choice and calculate their difference (Spectrum1-Spectrum2)
###### Input needed:
AveragedData.txt, AveragedReferenceData.txt of the same size, in textboxes: 'spectroscopy type'
###### Output:
FilenameOfYourAveragedData_Difference.txt

#### show points:
Import a data set with wavenumber/Raman shift/wavelength values (which you may create e. g. via the button 'point picking')
and highlight these positions in your latest matplotlib window as red vertical lines.
###### Input needed:
AnyData.txt, AnyData_points.txt
###### Output:
None

#### show graph:
Plots a data set of your choice.
###### Input needed:
AnyData.txt, in textboxes: 'spectroscopy type'
###### Output:
None

#### course plot:
Plots the temporal course of bands of your choice from your time-resolved data in 2D I vs t diagram.
###### Input needed:
RawData.txt (for DRIFTS: time values as a file with the same name as the raw data with the addition of _t ('FilenameOfYourRawData_t.txt'))
and AnyData_points.txt, in textboxes: 'spectroscopy type', 'number of spectra per period', 'number of periods' (choose 1 for averaged spectra),
'number of periods to cut off', for Raman: 'spike removing'
###### Output:
Courses.txt
  

