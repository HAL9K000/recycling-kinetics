# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 07:15:26 2020

@author: Koustav
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pan
import seaborn as sea
from scipy.optimize import curve_fit
import os


def f(x,a,b,c):
    return a+ b*np.exp(-c*x)


def plotter():
    
    os.chdir("Plot\MultiPlot")
    
    string="30 Hz"
    
    data=np.genfromtxt('Binder %s.csv' %(string), delimiter=",", skip_header=1)
    
    masterbinder={} #Empty dictionary to store all the values
    
    labels=["Trial Num","Time","Total Fluorescence","RRP Size","Endoctosed Pool Size","Fused Pool Size","Type"]
    
    for x in range(0, len(labels)-1):
        masterbinder[labels[x]]=data[:,x]
    #Created masterbinder.
    
    masterbinder['Type']=[]
    
    for x in range(0,len(data[:,6])):
        if data[x,6]== 0:
            masterbinder['Type'].append("Theoretical")
        else:
            masterbinder['Type'].append("FM143")
    
    jailkeeper= pan.DataFrame(masterbinder)
        
    g= sea.lineplot(x="Time", y="Total Fluorescence",hue="Type", estimator='mean', ci='sd' ,data=jailkeeper)
    plt.ylim(0,30)
        
    plt.savefig("%s Cumulative Plot.png" %(string), dpi=300)
        #plt.show()
        #plt.close()
        
    a=[]; fmean=[];fsd=[]
    for x in range(1,800):
        for y in range(0,len(masterbinder["Fused Pool Size"])):
            if masterbinder["Time"][y]==x:
                a.append(masterbinder["Total Fluorescence"][y])
        fmean.append(np.mean(a, dtype=np.float64))
        fsd.append(np.std(a, dtype=np.float64))
        
    print(fmean)
    print(fsd)
    
    #For 1 Hz     
    #popt, pcov = curve_fit(f, list(range(1,800)), fmean, sigma= fsd, bounds=([10.0,4.0,0.001],[30.0,20.0, 0.01]))
    
    #For 3 Hz
    #popt, pcov = curve_fit(f, list(range(1,800)), fmean, bounds=([10.0,11.5,0.006],[20.0,16.0, 0.015]))
    
    #For 10 Hz
    popt, pcov = curve_fit(f, list(range(1,800)), fmean, sigma= fsd, bounds=([5.0,15.0,0.01],[25.0,25.0, 0.03]))
    
    #For 30 Hz
    #popt, pcov = curve_fit(f, list(range(1,800)), fmean, sigma= fsd, bounds=([8.0,0,0.08],[15.0,25.0, 0.08]))
    
    xdata=np.array(list(range(0,800)))
        
    plt.plot(xdata, f(xdata, *popt), 'm--', label='Th Fit: %5.4f + %5.4f*e^(-%5.4fx)' % tuple(popt))
    #plt.savefig("%s Cumulative Plot Alt.png" %(str), dpi=300)
    plt.savefig("%s Best Fit Cumulative Plot.png" %(string), dpi=300)
    plt.show()
    plt.close()
    
    
plotter()
    
        
    