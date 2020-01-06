# -*- coding: utf-8 -*
"""
Created on Mon May 13 09:15:18 2019 by Jakob Weyel, TU Darmstadt

Als Schmankerl kann diese Version auch Fourierreihen in der Funktion
PSD_rechnen mit den Signalen falten. Dafür nur den oberen Wert von k bei
"#Kleinkram definieren" beliebig erhöhen!

Immer aufpassen, dass die Zeilen zum einschmuggeln eines periodischen Signals
in PSD_rechnen auskommentiert sind, falls nicht explizit erwünscht!

Richtig nervig: Nach Bandenauswahl und Banden_anzeigen muss die GUI noch per
destroy-Befehl geschlossen werden (oftmals auch per Hand), da sonst nichts im
Bildchen angezeigt wird. Irgendwie stören sich da tkinter und die matplotlib!

@author: jakob.weyel@tu-darmstadt.de
"""
import time

from tkinter import *
from tkinter import messagebox
from tkinter import filedialog

import numpy as np
from scipy import integrate as igr
import matplotlib.pyplot as plt
import pandas as pd

'''
_______________________________________________________________________________
ABFRAGE- (YES/NO) UND SPEICHERBUTTON
_______________________________________________________________________________
'''

def yesno(name, output, text): #Hier wird entschieden, ob die TRS-Spektren, PSD-Spektren oder die Bandenlagen mit in-Phase-Winkeln gespeichert werden
    
    msgbox = messagebox.askquestion ('Obacht!', text, icon = 'warning')
    if msgbox == 'yes':
       output.to_csv(name, sep = '\t', index = False)
       messagebox.showinfo('Ja, Mann!', 'gespeichert als: ' + name)
    else:
        messagebox.showinfo('Pech.','Nichts wurde gespeichert.')

'''
_______________________________________________________________________________
DATEIEN ÖFFNEN
_______________________________________________________________________________
'''

def myFileOpen(text): #Lädt die gewünschte Datei ein
	filename = filedialog.askopenfilename(title = text)
	return filename

'''
_______________________________________________________________________________
RECHENSPÄßE
_______________________________________________________________________________
'''

def TRS_rechnen():
    #Anzahl an Spektren pro Periode und die Dateipfade werden über die GUI eingelesen
    text = 'Dateipfad der Kat.-Spektren'
    name_mitgas = myFileOpen(text)
    
    text = 'Dateipfad der Ref.-Spektren'
    name_nurgas = myFileOpen(text)
    
    #Daten von Katalysator mit Gasphase und referenz mit Gasphase einlesen
    mitgas = np.genfromtxt(r''+name_mitgas, delimiter="\t")
    nurgas = np.genfromtxt(r''+name_nurgas, delimiter="\t")
    
    
    
    
    
    #Differenzspektren bilden, um Gasphase abzuziehen
    data = mitgas-nurgas
    data[:,0] = mitgas[:,0]
    
    #Zeitaufgelöste Spektren plotten
    plt.figure()
    plt.xlabel(r'$\tilde{\nu}$ / cm$^{-1}$')
    plt.ylabel('Int.')
    
    plt.plot(data[:,0],data[:,1:], linewidth = 0.1)
    
    plt.ylim(-np.amax(data[:,1:]), np.amax(data[:,1:]))
    
    #plot immer tight in page
    plt.tight_layout()
    
    #TRS-Spektren Speichern?
    TRS = ['TRS ' + str(i) for i in range(1,len(data[0,:]))]
    
    header = ['Wavenumber'] + TRS
    output = pd.DataFrame(data = data[:,:], columns = header)
    
    #Soll der Quatsch gespeichert werden?
    
    text = 'Sollen die TRS-Spektren als .txt gespeichert werden?' #Speicherabfrage
    name = name_mitgas.split('.') #Wird beim Speichern als Dateiname verwendet
    name = name[0] + '_TRS_Spektrum.txt'
    yesno(name, output, text)

    return
    
def PSD_rechnen(): #Macht aus trs- PSD-Spektren
    #Anzahl an Spektren pro Periode, rauszuschneidende Perioden, Phasenauflösung und die Dateipfade werden über die GUI eingelesen
    start = time.time()
    
    n_sp = int(myEntry01.get()) #Anzahl Spektren pro Periode
    raus_per = int(myEntry02.get()) #Rauszuschneidende Perioden
    raus_sp = n_sp*raus_per #Rauszuschneidende Spektren
    dphi = int(myEntry03.get()) #Phaseninkrement
    
    text = 'Dateipfad der Kat.-Spektren'
    name_mitgas = myFileOpen(text)
    
    name_t = name_mitgas.split('.')
    name_t = name_t[0] + '_t.dpt'
    
    text = 'Dateipfad der Ref.-Spektren'    
    name_nurgas = myFileOpen(text)
    
    #Daten von Katalysator mit Gasphase und KBr mit Gasphase einlesen
    mitgas = np.genfromtxt(r''+name_mitgas, delimiter="\t")
    t_inp = np.genfromtxt(r''+name_t, delimiter="\t")
    
    if name_nurgas!= '':
        nurgas = np.genfromtxt(r''+name_nurgas, delimiter="\t")
        
        #Differenzspektren bilden, um Gasphase abzuziehen
        data = mitgas-nurgas
        data[:,0] = mitgas[:,0]
        
    else :
        data = mitgas
    
    if raus_sp != 0:
        #Zeug abschneiden, was halt in Einpendelphase liegt
        data = np.delete(data, np.s_[1:raus_sp+1], axis = 1)
        t_inp = np.delete(t_inp,np.s_[-(raus_sp):],axis = 0)
        
    #Kleinkram definieren
    k = 1
    omega = 2*np.pi/t_inp[int(n_sp),0] #omega als 2pi/Periodendauer, Periodendauer aus Anzahl der Messungen bis zum zweiten GP-Wechsel
    phi = np.arange(0,361,dphi) #Phi ist die Phasenverschiebung, die bei Antwort des Systems auf GP-Wechsel erwartet wird
    
#    #Zum testen was periodisches einschmuggeln (nur für Notfälle!)
#    for l in range(0, len(data[0,1:]-1)):
#        data[10,1:] = (omega*t_inp[:,0]/1000000)
    
    Spektrum = np.zeros((len(data[:,0]),len(phi)+1))
    Spektrum[:,0] = mitgas[:,0]
    sos = Spektrum
    
    #Für alle vorgegebenen phi und Wellenzahlen die Fouriertrafo von Zeit -> Phase durchführen
    for k in np.arange(1,2): #Wenn man k höher wählt wird per Fourierreihe eine Rechteckfunktion reingewurschtelt
        for i in range(1,len(phi)+1):
            for j in range(0,len(data[:,0])):
                sos[j,i] = igr.trapz(data[j,1:]*(1/(2*k))*np.sin((2*k-1)*omega*t_inp[:,0]+phi[i-1]*2*np.pi/360))
    Spektrum[:,1:] = Spektrum[:,1:]+sos[:,1:]
    
    #Spektren plotten
    plt.figure()
#    for i in range(1,len(phi)+1):
    
    plt.plot(Spektrum[:,0],Spektrum[:,1:])#,label = str(phi[i-1]))
        
    #PSD-Daten in hübschen Dataframe verpacken
    phi_str = [str(item) + '°' for item in phi]
        
    header = ['Wavenumber'] + phi_str
    output = pd.DataFrame(data = Spektrum[:,:], columns = header)
    
    plt.ylim(-np.amax(Spektrum[:,1:]), np.amax(Spektrum[:,1:]))
    
    plt.xlabel(r'$\tilde{\nu}$ / cm$^{-1}$')
    plt.ylabel(r'-lg($R$)')
    
    #Legende rein
    plt.legend(phi, title = r'$\phi^\mathrm{PSD}$ / °', loc = 'upper right')
    
    end = time.time()
    
    #Soll der Quatsch gespeichert werden?
    
    text = 'Sollen die PSD-Spektren als .txt gespeichert werden?' #Speicherabfrage
    name = name_mitgas.split('.') #Wird beim Speichern als Dateiname verwendet
    name = name[0] + '_PSD_Spektrum_' + str(raus_per) + '_Perioden_raus_' + str(dphi) + '_dphi.txt'
    yesno(name, output, text)
    print('Runtime PSD: ' + str(end-start) +' s.')

    return

def Bandenauswahl():
    fig = plt.figure()
    text = 'Dateipfad der PSD-Spektren'
    name_psd = myFileOpen(text)
    psd = np.genfromtxt(r''+name_psd, delimiter="\t", skip_header=1)
    #Alle PSD-Spektren gleichzeitig plotten
    plt.clf()
    plt.xlabel(r'$\tilde{\nu}$ / cm$^{-1}$')
    
    for i in range(1, psd.shape[1]-1):
        plt.plot(psd[0:,0],psd[0:,i], 'o', picker = 3)
#        plt.xlabel(r'$\tilde{\nu}$ / cm$^{-1}$')
        plt.ylim(-np.amax(psd[:,1:]), np.amax(psd[:,1:]))
        plt.show()
        
        bands = []
        
        def onpick(event): #Das Moped erlaubt es, die aktuelle Position des Mauszeigers bei anklicken von Datenpunkten in array zu speichern
            thisline = event.artist
            xdata = thisline.get_xdata()
            ind = event.ind
            points = xdata[ind]
            bands.append(points[int(len(points)/2)]) #Array von Bandenpositionen wird um aktuellen Klick verlängert
            
            #Hier wird gespeichert, sobald ein Datenpunkt angeklickt wird!!!
            bands_new = sorted(set(bands)) #sorted sortiert aufsteigend, set soriert und entfernt doppelte
            name = name_psd.split('.')
            bands_new = np.array(bands_new)
            np.savetxt(name[0] + '_Banden.txt',bands_new, delimiter = '\t')
            
    fig.canvas.mpl_connect('pick_event', onpick)
    myPSDgui.destroy() #Irgendwas an tkinter stört die matplotlib, deswegen wird die GUI zugemacht (Zur Not alles außer Plotfenster per Hand schließen!)
    return

def Banden_anzeigen():
    text = 'Dateipfad der Bandenlagen'
    name_Banden = myFileOpen(text)
    Banden = np.genfromtxt(r''+name_Banden, delimiter="\n")
    for ii in np.arange(0,len(Banden)):
        ymin,ymax = plt.gca().get_ylim()
        plt.plot([Banden[ii],Banden[ii]],[ymin,ymax],'r', linewidth = 0.5)
    myPSDgui.destroy() #Irgendwas an tkinter stört die matplotlib, deswegen wird die GUI zugemacht (Zur Not alles außer Plotfenster per Hand schließen!)
    return

def in_Phase_Winkel():
    #Erstmal die PSD-Spektren und die Bandenlagen einlesen
    text = 'Dateipfad der PSD-Spektren'
    name_psd = myFileOpen(text)    
    psd_Spektren = pd.read_csv(r''+name_psd, delimiter="\t")
    #Hier muss gerundet werden, da pandas mehr Nachkommastellen einlädt als numpy
    psd_Spektren.Wavenumber = pd.Series([round(val, 5) for val in psd_Spektren.Wavenumber], index = psd_Spektren.index)
    
    text = 'Dateipfad der Bandenlagen'
    name_Banden = myFileOpen(text)
    Bandenlagen = np.genfromtxt(r''+name_Banden, delimiter="\n")
    Bandenlagen = Bandenlagen[::-1] #Muss invertiert werden, sonst passiert Unfug bei der Zuordnung!
    
    '''Hier werden die vorgegebenen Bandenlagen mit den Wellenzahlen im
    Spektrum verglichen und der nächstgelegene Wert weiterverwendet'''
    i = 0
    for val in Bandenlagen:
        Bandenlagen[i] = min(psd_Spektren.Wavenumber, key=lambda x:abs(x-val))
        i = i+1
        
    #Zeiten einlesen, um später von Phasenwinkel auf die Zeit bis Signalmaximum umzurechnen    
    
    text = 'Dateipfad der Zeitwerte'
    name_t = myFileOpen(text)
    t_inp = np.genfromtxt(r''+name_t, delimiter="\t")
    n_sp = int(myEntry01.get())
    tges = t_inp[n_sp-1,0]
    
    #Die Zeilen, die zu gewählten Wellenzahlen gehören (mit Bande) vom restlichen Unfug trennen!
    phi_bei_Banden = psd_Spektren[psd_Spektren.Wavenumber.isin(Bandenlagen)]
    phi_bei_Banden = phi_bei_Banden.iloc[:,1:] #Blöde wellenzahlspalte loswerden
    
    #Phasenwinkel auslesen, bei dem das zur jeweiligen Wellenzahl gehörige Maximum liegt.
    Wmax = phi_bei_Banden.idxmax(axis = 1)
    Wmax = np.array(Wmax.values,dtype = str)
    Wmax = np.core.defchararray.replace(Wmax,'°','')
    Wmax = np.array(Wmax,dtype = int)
    
    #Aus Maximumwinkel die Zeit ausrechnen!
    tmax = (360-Wmax)/360*tges
    
    #Outputarrays der Wellenzahl und Zeiten werden auf erste Nachkommastelle gerundet
    Bandenlagen = np.around(Bandenlagen,1)
    tmax = np.around(tmax,1)
    
    output = pd.DataFrame({'Wavenumber': Bandenlagen, 'Phi_max / °': Wmax, 't / s mit t_Per. = ' + str(tges) + ' s': tmax})
    
    #Hier wird gespeichert!
    text = 'Sollen die Bandenlagen mit den in-Phase-Winkeln als .txt gespeichert werden?' #Speicherabfrage
    name = name_psd.split('.') #Wird beim Speichern als Dateiname verwendet
    name = name[0] + '_Banden_iPW.txt'
    yesno(name, output, text)
  
def Bilder_anzeigen():
    text = 'Dateipfad der PSD-Spektren'
    name_psd = myFileOpen(text)
    psd = np.genfromtxt(r''+name_psd, delimiter="\t", skip_header=1)
    
    #Alle PSD-Spektren gleichzeitig plotten
    plt.figure()#figsize=(5,3), dpi=200)
    
    plt.plot(psd[:,0],psd[:,1:])
    
    plt.xlabel(r'$\tilde{\nu}$ / cm$^{-1}$')
    plt.ylim(-np.amax(psd[:,1:]), np.amax(psd[:,1:]))
    
    # Immer aufpassen, ob Legende und Achsenbeschriftuungen Sinn machen!
    phi = np.arange(0,361,30)
    plt.legend(phi, title = '$\phi^\mathrm{PSD}$ / °', loc = 'upper right')

#    plt.yticks([],[])
    plt.ylabel('-lg($R$)')

    return

def Spur():
    #Spektrensätze, Zeit und Bandenlagen einladen
    n_sp = int(myEntry01.get()) #Anzahl Spektren pro Periode
    raus_per = int(myEntry02.get()) #Rauszuschneidende Perioden
    raus_sp = n_sp*raus_per #Rauszuschneidende Spektren
    
    text = 'Dateipfad der Kat.-Spektren'
    name_mitgas = myFileOpen(text)
    
    name_t = name_mitgas.split('.')
    name_t = name_t[0] + '_t.dpt'
    
    text = 'Dateipfad der Ref.-Spektren'
    name_nurgas = myFileOpen(text)
    
    text = text = 'Dateipfad der Bandenlagen'
    name_Banden = myFileOpen(text)        
    
    #Daten von Katalysator mit Gasphase, Referenz mit Gasphase und Bandenlagen einlesen
    mitgas = np.genfromtxt(r''+name_mitgas, delimiter="\t")
    t_inp = np.genfromtxt(r''+name_t, delimiter="\t")
    Banden = np.genfromtxt(r''+name_Banden, delimiter="\n")
    
    if name_nurgas!= '':
        nurgas = np.genfromtxt(r''+name_nurgas, delimiter="\t")
        
        #Differenzspektren bilden, um Gasphase abzuziehen
        data = mitgas-nurgas
        data[:,0] = mitgas[:,0]
        
    else :
        data = mitgas
    
    if raus_sp != 0:
        #Zeug abschneiden, was halt in Einpendelphase liegt
        data = np.delete(data, np.s_[1:raus_sp+1], axis = 1)
        t_inp = np.delete(t_inp,np.s_[-(raus_sp):],axis = 0)
    
    t_1Spektrum = t_inp[1,0]-t_inp[0,0] # Zeit, die die Aufnahme eines Spektrums braucht als Differenz zweier konkreter Zeitwerte ausrechnen
    print('Ein Spektrum brauchte ' + str(t_1Spektrum) + ' s')
    
    i = 0
    for val in Banden:
        Banden[i] = min(data[:,0], key=lambda x:abs(x-val)) # Wenn Wellenzahlen aus "Banden" und data[:,0] nicht 100pro übereinstimmen, wird der nächstliegende Wert mit der geringsten Abweichung genommen.
        i = i+1
        
    #Jetzt wird die Spur an den jeweiligen Bandenlagen geplotted
    pos = np.zeros(len(Banden))
    for i in np.arange(0,len(Banden)):
        if i%10 == 0:# and i != 0: #Macht alle zehn Linien ein neues Plotfenster auf, da die Farben sich sonst wiederholen
            #Phasen zur Unterscheidung einfärben (in allen außer dem letzten Bildchen)
            plt.figure()#figsize=(10,6), dpi=80)
            for j in np.arange(min(t_inp[:,0]),max(t_inp[:,0]/60),n_sp*t_1Spektrum/60): #durch 60 teilen, um von s zu min zu kommen
                plt.axvspan(j, j+n_sp/2*t_1Spektrum/60, facecolor='k', alpha=0.25)
        
        sos = np.where(data[:,0] == Banden[i])
        pos[i] = sos[0]
        
        plt.plot(t_inp[:,0]/60,data[int(pos[i]),1:], label = str(np.around(Banden[i],1)))
        plt.legend(title = r'$\tilde{\nu}$ / cm$^{-1}$')
        plt.xlabel('$t$ / min')
        plt.ylabel('-lg($R$)')
    
'''
_______________________________________________________________________________
GUI-Gedöns
_______________________________________________________________________________
'''

myPSDgui = Tk()
myEntry01 = StringVar()
myEntry02 = StringVar()
myEntry03 = StringVar()

#myPSDgui.geometry('400x200+200+100')

myPSDgui.title('Hier kommen PSD-Spektren raus!')

myLabel01 = Label(myPSDgui, text = 'Wie viele Spektren wurden pro Periode aufgenommen (ox. + red.)?').pack()

myEntry01 = Entry(myPSDgui, textvariable = myEntry01)
myEntry01.insert(END,'40')
myEntry01.pack() #Bei Entryfeld muss .pack() in neue Zeile, sonst ist der Output = None!

myLabel02 = Label(myPSDgui, text = 'Wie viele Perioden als Einpendelphase abschneiden?').pack()

myEntry02 = Entry(myPSDgui, textvariable = myEntry02)
myEntry02.insert(END,'0')
myEntry02.pack() #Bei Entryfeld muss .pack() in neue Zeile, sonst ist der Output = None!

myLabel03 = Label(myPSDgui, text = 'Wie fein soll die Phasenauflösung sein?').pack()

myEntry03 = Entry(myPSDgui, textvariable = myEntry03)
myEntry03.insert(END,'30')
myEntry03.pack() #Bei Entryfeld muss .pack() in neue Zeile, sonst ist der Output = None!


myButton01 = Button(myPSDgui, text = 'TRS -> PSD', command = PSD_rechnen).pack()

myLabel04 = Label(myPSDgui, text = 'Wenn kein Bildchen zu sehen ist, alle Fenster von tkinter schließen! \n Sobald ein Datenpunkt angeklickt wird, werden alle bisher \n angeklickten Datenpunkte in eine Datei geschrieben:').pack()

myButton02 = Button(myPSDgui, text = 'Bandenauswahl', command = Bandenauswahl).pack()

myLabel05 = Label(myPSDgui, text = 'Hier wird der Satz an TRS-Spektren bearbeitet:').pack()

myButton03 = Button(myPSDgui, text = 'Referenzabzug', command = TRS_rechnen).pack()

myLabel06 = Label(myPSDgui, text = 'Bandenlagen in aktuell ausgewähltem Bildchen anzeigen:').pack()

myButton04 = Button(myPSDgui, text = 'Banden anzeigen', command = Banden_anzeigen).pack()

myLabel07 = Label(myPSDgui, text = 'Winkel für in-Phase-Spektren bestimmen:').pack()

myButton05 = Button(myPSDgui, text = 'in-Phase-Winkel', command = in_Phase_Winkel).pack()

myLabel08 = Label(myPSDgui, text = 'Einfach nur Bilder gucken (PSD oder TRS oder sonstwas):').pack()

myButton06 = Button(myPSDgui, text = 'Spektrensatz anzeigen', command = Bilder_anzeigen).pack()

myLabel09 = Label(myPSDgui, text = 'Spurplots ausgewählter Banden anzeigen:').pack()

myButton07 = Button(myPSDgui, text = 'Spurplot', command = Spur).pack()


myPSDgui.mainloop()