# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 16:57:42 2020

@author: Koustav
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import os
import random as ran
import pandas as pan
import seaborn as sea
from scipy.optimize import curve_fit
from General_Mark_II import GeneralTwo

'''
Used to make multi-plots using seaborn for each trophic level

'''


'''

For 1 Hz:
    
    self.x=5
    self.x2= 0.05
    self.lambd= -0.0056
    self.alph=-0.0019
    self.alphpr= -0.0025
    
    Retain Eq:
        0.66587+ 0.2391*math.exp(self.ep*self.t)
        
        (Derived From SVFR)
        
For 3 Hz:
    
    self.x=5
    self.x2=0.05
    self.lambd=-0.0132
    self.alph= -0.0031
    self.alphr= -0.0025
    self.ppr=0.5
    
    By Cubic Method:
        
       self.alph= -0.00375 OR -0.02719
       
    
    Retain Eq:
        0.5823+ 0.4177*math.exp(self.ep*self.t)
        (Derived From SVFR)
        
    What actually works:
        0.65+ 0.35*math.exp(self.ep*self.t)         (Inexplicable)
        
        self.x2=0.12
        
For 10 Hz:
    
    self.x=5
    self.lambd=-0.0310
    self.alph= -0.0114
    self.alphr= -0.0025
    
    Using Cubic Eq:
    
        self.alph= -0.02809
    
    Retain Eq:
        0.3252+ 0.6748*math.exp(self.ep*self.t)
        (Derived From SVFR)
        
For 30 Hz:
    
    self.x=10
    
    self.lambd= -0.0621
    self.alph= -0.0401
    self.alphr= -0.0025
    
    Retain Eq:
        0.49079+ 0.50920*math.exp(self.ep*self.t)
        (Derived From SVFR)
        
    Using Cubic Eq:
    
        self.alph= -0.05549
    
        
    
'''



class MultiPlotter(GeneralTwo):
    
    
    def __init__(self):
        
        self.xrange=np.array([5,6,7,8,9,10])  #Possible values that RRP size can attain
        self.x=ran.choice(self.xrange)        #Choose one RRP size at random from xrange
        self.x=10
        self.trps=30                        #Size of TRP
        
        
        
        
        self.ratepicker()       #Input a frequency and have all the relevant parameters set for you.
        
        'Function above called in from the class GeneralTwo'
        
        
        self.tracker=np.ones(self.trps)                     #Setting up tracker for individual vesicles with each element containing initial fluorescence.
        self.t=0              #Initial time
        self.dt=1          # Chosen as time step.
        self.pref=np.zeros([self.trps,4])
        '''Keeps track of changes taking place in vesicles, with 1st column storing kind of fusion 
        ["100" for k&r, "200" for full fusion, "300" for endocytosed ("301" for K&R, "302" for full fusion), 
        "400" for RRP, "500" for Active RCP, "0" for RCP unfused], 
        2nd column noting exact time of this occuring, the 3rd one noting the number of such consecutive fusions
        and the 4th storing a random number that determines the probability of transfer of each vesicle from one to another pool'''
        self.pref[0:self.x,0]=400
        self.eptrack=self.trps-self.x                    #Notes number of vesicles in RCP at the beginning
        self.fsum=[]                        #Measures total fluorescence at each time step.
        self.rrpsiz=[]
        self.endsiz=[]
        self.fussiz=[]
        self.rcpsiz=[]
        self.ftime=[]
        
        
        self.fmean=[] #Stores fluorescent means for all fifty iterations
        self.fsd=[]     #
        
        self.masterbinder={} #Empty dictionary that stores a lot of collated data
        
        
        for r in range(0,self.x):
            self.pref[r,3]= ran.random() #Generating random number for these beauts.
            
        self.sanitizer()
        
        
    def masterdeck(self):
        
        self.counter=0
        
        for x in range(0,35):
            self.counter+=1
            #50 simulations to be performed in total.
            self.controlpanel()
            print "Hogwarts"
            self.counter+=1
            
        self.statistics()  #Makes requisite plots at the end.
        
            
            
    def accountancy(self):   #To create a masterbinder that stores all the relevant data generated in each trial.
    
    
        if "Trial Num" not in self.masterbinder:
            self.masterbinder["Trial Num"]=[]
            self.masterbinder["Time"]=[]
            self.masterbinder["Total Fluorescence"]=[]
            self.masterbinder["RRP Size"]=[]
            self.masterbinder['Endoctosed Pool Size']=[]
            self.masterbinder['Fused Pool Size']=[]
            self.masterbinder['Type']=[]
            #Creating these entries in master-binder for first time 
            
        self.masterbinder["Total Fluorescence"].extend(self.fsum)
        self.masterbinder["RRP Size"].extend(self.rrpsiz)
        self.masterbinder["Endoctosed Pool Size"].extend(self.endsiz)
        self.masterbinder["Fused Pool Size"].extend(self.endsiz)
        self.masterbinder["Time"].extend(list(range(0,self.t)))
        
        for x in range(0,self.t):
            self.masterbinder['Type'].append(0)
            self.masterbinder['Trial Num'].append(self.counter)
            
        
        
        self.sanitizer()
        #Sanitise at the end of every trial run
        
    def statistics(self):
        #Crunches the numbers, makes all those fancy plots.
            
        os.chdir("Plot")
        if( os.path.isdir(self.str)==False):
            os.mkdir(self.str)
        os.chdir(self.str)
        fm143=np.genfromtxt("%s FM143.csv" %(self.str), delimiter=',')
        fm143[:,1]*=30
        
        
        a=[]
        for x in range(0,800):
            for y in range(0,len(self.masterbinder["Fused Pool Size"])):
                if self.masterbinder["Time"][y]==x:
                    a.append(self.masterbinder["Total Fluorescence"][y])
            self.fmean.append(np.mean(a, dtype=np.float64))
            self.fsd.append(np.std(a, dtype=np.float64))
        
        
        a,b =fm143.shape
        
        #Stores data related to FM143 dye.
        self.masterbinder["Time"].extend(fm143[:,0])
        self.masterbinder["Total Fluorescence"].extend(fm143[:,1])
        
        
        
        
        for x in range(0,a):
            self.masterbinder['Type'].append(1)         # 1 represents type FM143
            self.masterbinder["RRP Size"].append(-1)
            self.masterbinder["Endoctosed Pool Size"].append(-1)
            self.masterbinder["Fused Pool Size"].append(-1)
            self.masterbinder['Trial Num'].append(-1)
            
        os.chdir("../")
        if( os.path.isdir("MultiPlot")==False):
            os.mkdir("MultiPlot")
        os.chdir("MultiPlot")
        
        a=len(self.masterbinder["Total Fluorescence"])
        
        datasaver=np.zeros([a,7])
        
        datasaver[:,0]=self.masterbinder['Trial Num']
        datasaver[:,1]=self.masterbinder['Time']
        datasaver[:,2]=self.masterbinder['Total Fluorescence']
        datasaver[:,3]=self.masterbinder['RRP Size']
        datasaver[:,4]=self.masterbinder['Endoctosed Pool Size']
        datasaver[:,5]=self.masterbinder['Fused Pool Size']
        datasaver[:,6]=self.masterbinder['Type']
        
        header_text = "Trial Num,Time,Total Fluorescence,RRP Size,Endoctosed Pool Size,Fused Pool Size,Type"
        
        np.savetxt('Binder %s.csv' %(self.str), datasaver, delimiter=",", header=header_text,)
        
        '''jailkeeper= pan.DataFrame(self.masterbinder)
        
        g= sea.lineplot(x="Time", y="Total Fluorescence",hue="Type", estimator='mean', ci='sd' ,data=jailkeeper)
        plt.ylim(0,30)
        
        plt.savefig("%s Cumulative Plot.png" %(self.str), dpi=300)
        #plt.show()
        #plt.close()'''
        
        
        
        popt, pcov = curve_fit(self.f, list(range(0,800)), self.fmean, sigma= self.fsd, bounds=([5.0,0,0],[25.0,25.0, 0.05]))
        xdata=np.array(list(range(0,800)))
        
        plt.plot(xdata, self.f(xdata, *popt), 'm--', label='Th Fit: %5.4f + %5.4f*e^(-%5.4fx)' % tuple(popt))
        plt.savefig("%s Cumulative Plot Alt.png" %(self.str), dpi=300)
        plt.show()
        plt.close()
        
        
        
        
        
        
        
        
        
        
        
        
            
    def sanitizer(self): 
            # Santitizes and resets all variables to their default (initial) values. Should only be used at the beginning and ends of runs.
            
            self.x2=ran.choice(self.x2range())            #Chosen extent of full scale fusion
            #self.x2=0
        
            self.x1= 1 - self.x2                      #Chosen extent of kiss and run
            
            self.tracker=np.ones(self.trps)                     #Setting up tracker for individual vesicles with each element containing initial fluorescence.
            self.t=0              #Initial time
            
            self.pref=np.zeros([self.trps,4])
            
            self.pref[0:self.x,0]=400
            self.eptrack=self.trps-self.x                    #Notes number of vesicles in RCP at the beginning
            self.fsum=[]                        #Measures total fluorescence at each time step.
            self.rrpsiz=[]
            self.endsiz=[]
            self.fussiz=[]
            self.rcpsiz=[]
            self.ftime=[]
            
            for r in range(0,self.x):
                self.pref[r,3]= ran.random() #Generating random number for these beauts.
                
                
obj=MultiPlotter()    
if __name__ == "__main__":
    
    obj.masterdeck()
            
            