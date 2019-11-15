# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 01:35:27 2019

@author: Koustav
"""

from General_Mark_I.py import GeneralOne
import numpy as np
import math
import matplotlib.pyplot as plt
import os
import random as ran
from scipy.optimize import curve_fit

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
    
        self.alph= -0.05549 '''
    
        

class Generale(GeneralOne):
    
    def __init__(self):
        
        self.xrange=np.array([5,6,7,8,9,10])  #Possible values that RRP size can attain
        self.x=ran.choice(self.xrange)        #Choose one RRP size at random from xrange
        self.x=5
        self.trps=30                        #Size of TRP
        self.x2_low=0.10                    #Lower limit on full scale fusion (but why such a value??)
        self.x2_stepsize=0.01
        self.x2_high=0.20                   #Upper limit on full scale fusion
        self.x2=ran.choice(GeneralOne.x2range(self))            #Chosen extent of full scale fusion
        #self.x2=0.05
        print(self.x2)
        #self.x2= 0.15
        self.x1= 1 - self.x2                      #Chosen extent of kiss and run
        self.gamma= -2.0                                    #Look up the model for meanings of specific parameters
        self.delta=-1.0/120                             #Stevens & Murphy, 2018
        #self.deltap=
        self.lambd= -0.0310
        self.alph=-0.02809
        self.alphpr= -0.0025
        self.beta= -0.1667
        self.ppr=0.5
        self.epold= ((self.gamma+2*self.delta)*self.lambd)/(self.gamma+self.lambd)
        print "Old Epsilon:\t%f" %(self.epold)
        self.ep= (self.lambd/(self.gamma*(self.gamma+ self.lambd)))*((self.gamma)**2 + (self.ppr*self.delta)*(self.gamma+self.lambd))
        print "Epsilon:\t%f" %(self.ep)
        self.deprt= -math.log(2)/2.5                        #Chosen rate of departitioning of fm1-43 dye ( based on 2.5 s halftime)
        self.tracker=np.ones(self.trps)                     #Setting up tracker for individual vesicles with each element containing initial fluorescence.
        self.t=0              #Initial time
        self.dt=1          # Chosen as time step.
        self.pref=np.zeros([self.trps,4])
        '''Keeps track of changes taking place in vesicles, with 1st column storing kind of fusion 
        ("100" for k&r, "200" for full fusion, "300" for endocytosed, "400" for RRP, "500" for Active RCP, "0" for RCP unfused), 
        2nd column noting exact time of this occuring and the 3rd one noting the number of such consecutive fusions.'''
        self.pref[0:self.x,0]=400
        self.eptrack=self.trps-self.x                    #Notes number of vesicles in RCP at the beginning
        self.fsum=[]                        #Measures total fluorescence at each time step.
        self.rrpsiz=[]
        self.endsiz=[]
        self.ftime=[]
        self.str="10 Hz"             #For the purposes of making documentation easier.
        
