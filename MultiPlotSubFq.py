# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 21:14:22 2020

@author: Koustav
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import os
import random as ran
import pandas as pan
import seaborn
from scipy.optimize import curve_fit
from Alt_Mark_0_3 import PointthreeHz


''' Used to generate and output multiple trial data for 0.3 Hz sub-case '''

class MultiPlotterSubFq(PointthreeHz):
    
    
    def __init__(self):
        
        self.xrange=np.array([5,6,7,8,9,10])  #Possible values that RRP size can attain
        self.x=ran.choice(self.xrange)        #Choose one RRP size at random from xrange
        self.x =5
        self.trps=30                        #Size of TRP
        
        self.str="0.3 Hz"
        
        
        
        
        self.ratepicker()       #All the relevant parameters are set for you.
        
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
        
        
        self.fmean=[] #Stores fluorescent means for all  trials
        self.fsd=[]     #Stores fluorescent SDs for all  trials
        
        self.fm143sd= 0.5 #Check statistics() function. Used to define spread around experimental fm143 data.
        
        self.masterbinder={} #Empty dictionary that stores a lot of collated data
        
        self.trials = int(raw_input("Enter number of trials that you want:\t"))
        # Number of trials to be initiated
        
        for r in range(0,self.x):
            self.pref[r,3]= ran.random() #Generating random number for these beauts.
            
        for r in range(self.x,self.trps):
            self.pref[r,3]= ran.random() #Generating random number for these beauts
            
        self.sanitizer()
        
        
    def masterdeck(self):
        
        self.counter=0
        
        for x in range(0,self.trials):
            self.counter+=1
            #50 simulations to be performed in total.
            self.controlpanel()
            print "Hogwarts"
            
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
        
        for x in range(0,a):
            fm143_gauss = np.random.normal(loc=fm143[x,1], scale= self.fm143sd, size= self.trials)
            # fm143_gauss stores 35 samples drawn from a random distribution with mean = fm143[x,1] & sd = scale.
            print "Size of fm143_gauss is:\t %d" %(fm143_gauss.size)
            print fm143_gauss
            for y in range(0, fm143_gauss.size):
                self.masterbinder["Time"].append(fm143[x,0])
                #Stores time point associated with FM143 datapoint
            self.masterbinder["Total Fluorescence"].extend(fm143_gauss)
            #Stores Gaussian data point in it's full glory, alongside it's spread
        
        
        b= a*self.trials        #Used to fill in the remaining elements in the masterbinder for FM143 data.
        
        for x in range(0,b):
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
        
        np.savetxt('Binder De Novo %s.csv' %(self.str), datasaver, delimiter=",", header=header_text,)
        
        jailkeeper= pan.DataFrame(self.masterbinder)
        
        print seaborn.__version__
        
    def sanitizer(self): 
            # Santitizes and resets all variables to their default (initial) values. Should only be used at the beginning and ends of runs.
            
            ran.seed()
            
            self.x2=ran.choice(self.x2range())            #Chosen extent of full scale fusion
            self.x2=0       #NO FULL SCALE FUSION
        
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
                
            for r in range(self.x,self.trps):
                self.pref[r,3]= ran.random() #Generating random number for these beauts
                
obj=MultiPlotterSubFq()    
if __name__ == "__main__":
    
    obj.masterdeck()