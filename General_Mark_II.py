# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 15:32:29 2020

@author: Koustav
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import os
import random as ran
from scipy.optimize import curve_fit
from General_Mark_I import GeneralOne


'''
Major Difference with Version I is that we are now considering seperate rates (and a route) for Full Scale Fusion vesicle.

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


class GeneralTwo(GeneralOne):
    
    def __init__(self):
        
        self.xrange=np.array([5,6,7,8,9,10])  #Possible values that RRP size can attain
        self.x=ran.choice(self.xrange)        #Choose one RRP size at random from xrange
        self.x=10
        self.trps=30                        #Size of TRP
        
        
        
        
        self.ratepicker()       #Input a frequency and have all the relevant parameters set for you.
        
        
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
        
        for r in range(0,self.x):
            self.pref[r,3]= ran.random() #Generating random number for these beauts.
        
    def lrcp_rrp(self):                     #Determines how many vesicles have moved out from unfused RCP to RRP based on kinetic data.
    
        print "Mozart"

        retain= (self.trps-self.x)*(self.a+ self.b*math.exp(self.ep*self.t))
        #print "Retain:", retain
        # retain documents how many TRP vesicles haven't made the rcp-rrp transit at a particular time.
        if (self.eptrack-retain >=1):
            m=self.trps-self.eptrack         #Notes the vesicle number that will undergo transit
            k=int(self.eptrack-math.ceil(retain))        #In the event that there are more than one vesicles making the transit.
            print "k:",k
            for i in range(0,k):
                a=int(m+i)
                print "a:",a
                self.pref[a,0]=400          #Change tag of vesicle to RRP
                self.pref[a,1]=self.t       #Change timing to reflect entry to RRP\
                self.pref[a,3]=ran.random()
            self.eptrack=math.ceil(retain)       #Update eptrack to reflect current number of unfused vesicles
            for x in range(0,self.trps):
                print self.pref[x,0],
                
        for y in range(self.trps):
            if (self.pref[y,0]==0 and self.pref[y,1]>0):
                #Vesicle has been undocked to  Latent RCP.
                m=self.t-self.pref[y,1]         #Notes the time spent by vesicle in RCP.
                
                ch=self.a+ self.b*math.exp(self.ep*self.t)
                if (self.pref[y,3]>ch):
                    #Vesicle is primed and docked at RRP
                    self.pref[y,0]=400
                    self.pref[y,1]=self.t   #Updating time of entry to RRP.
                    self.pref[y,3]=ran.random()
        
        
    def fusion(self):                   #Goes into Narcos mode, checking for and updating vesicle parameters as they fuse with membrane.
        for y in range(self.trps):
            if (self.pref[y,0]==400):       #Checkes whether a vesicle is in the RRP.
                m=self.t - self.pref[y,1]        #Notes the amount of time a vesicle has spent in the RRP
                ch=math.exp((self.sdelta*m))
                #ch1=math.exp((self.sdeltap*m))
                if( self.pref[y,3] > ch):
                    #Vesicle is undocked from RRP......
                    self.pref[y,3]=ran.random()
                    if ( ran.random() <=0.5):
                        #...... To Latent RCP
                        self.pref[y,0]=0
                        self.pref[y,1]=self.t   #Updating time of entry to RCP.
                        print "CB'"
                        
                    else:
                        #.... To Active RCP.
                        self.pref[y,0]=500
                        self.pref[y,1]=self.t   #Updating time of entry to Active RCP.
                        print "CB"
     
                    continue #Move to next iteration
                    
                    
             
                    
                    
                    
                ch=math.exp((self.gamma*m))
                if( self.pref[y,3] > ch):              #Fusion has taken place.
                   print "CD!"
                   chec=ran.random()
                   if (chec <= self.x2):           #Full scale fusion has taken place.
                       print "Bo"
                       self.pref[y,0]=200
                       self.tracker[y]=0                #Setting fluorescence parameter to be 0
                       self.pref[y,3]=ran.uniform(0.3, 0.6)      
                       #Updating random-number stored so that the vesicle will undergo endocytosis between 15-35 s
                   elif (chec > self.x2):            # K&R fusion has taken place.
                        self.pref[y,0]=100
                        #self.tracker[y]*= math.exp(self.deprt*0.1*9)
                        self.tracker[y]*= math.exp(self.deprt*0.1*ran.choice(range(6,10)))
                        
                        self.pref[y,3]=ran.uniform(0.42, 1) 
                        #K&R vesicle  must not spend more than 2 second docked at the membrane
                        
                   self.pref[y,1]= self.t     #Update the time parameter to reflect fusion
                   
                   self.pref[y,2]+=1        #Notes the number of times fusion has taken place.
            else: continue
        
    def endo(self):   #Checks whether vesicle is endocytosed and updates parameters accordingly.
        for y in range(self.trps):
            m=self.t - self.pref[y,1]
            ch1=math.exp((self.delta*m))
            ch2=math.exp(((-math.log(2)/20.0)*m))
            if( self.pref[y,0]==100 and self.pref[y,3] > ch1):
                #K&R vesicle is endocytosed.
                self.pref[y,0]=301          
                self.pref[y,1]+=m         #Time of entering endocytosed compartment is updated.
                self.pref[y,3]=ran.random()
            elif (self.pref[y,0]==200 and self.pref[y,3] > ch2):
                #Full fusion vesicles are endocytosed
                self.pref[y,0]=302          
                self.pref[y,1]+=m         #Time of entering endocytosed compartment is updated.
                self.pref[y,3]=ran.random()
                
    def priming(self):  #Checks whether endocytosed vesicle is primed and updates parameters accordingly.
        print "YoloBobo"
        for y in range(self.trps):
            if (self.pref[y,0]==301):
                # K&R Vesicle is in endocytosed compartment.
                m=self.t-self.pref[y,1]         #Notes the time spent by vesicle as endocytosed.
                ch2=math.exp((self.alphpr*m))
                print "Alpha Prime 1%f" %(self.alph)
                print "Alt Alpha Prime %f" %(self.alp)
                #self.alph=self.alp
                ch1=math.exp((self.alph*m))
                
                
                if (ran.random()>ch1):
                    #Vesicle is primed and docked at ARCP
                    self.pref[y,0]=500
                    self.pref[y,1]=self.t   #Updating time of entry to ARCP.
                    #self.pref[y,3]=ran.random()
                    print "Blue Moon"
                    print "Alpha Prime 2%f" %(self.alph)
                    
                if (ran.random()>ch2):
                    #Vesicle is primed and docked at RRP
                    self.pref[y,0]=400
                    self.pref[y,1]=self.t   #Updating time of entry to RRP.
                    #self.pref[y,3]=ran.random()
                    continue
                
            if (self.pref[y,0]==302):
                print "Frank"
                #Full scale fusion vesicle is in endocytosed compartment.
                m=self.t-self.pref[y,1]         #Notes the time spent by vesicle as endocytosed.
                ch3=math.exp((self.alpha_full*m))
                if (ran.random()>ch3):
                    #Full scale fusion vesicle moves to ARCP
                    self.pref[y,0]=500
                    self.pref[y,1]=self.t   #Updating time of entry to ARCP.
                    #self.pref[y,3]=ran.random()
                    print "Frank Sinatra"
                    
    def acrp_rrp(self):  #Checks whether vesicle in Act RCP gets docked at the RRP
        for y in range(self.trps):
            if (self.pref[y,0]==500):
                #Vesicle is in Active RCP
                m=self.t- self.pref[y,1]
                ch=math.exp(self.beta*m)
                if (self.pref[y,3]> ch):
                    #Vesicle gets docked at RRP
                    self.pref[y,0]=400
                    self.pref[y,1]=self.t #Updating time of entry to RRP.
                    self.pref[y,3]=ran.random()
                    print "Holy Guacamole"
           
        
    def ratepicker(self):
        
        f= int(raw_input('Enter input frequency (choose b/w 1, 3, 10 or 30 Hz):\t'))
        
        self.x2_low=0.02                    #Lower limit on full scale fusion (but why such a value??)
        self.x2_stepsize=0.01
        self.x2_high=0.06                   #Upper limit on full scale fusion
        
        self.x=5        #Size of RRP
        self.trps=30    #Total size of TRP
        
        '''For a detailed understanding of the various rates, please look at the PDF documentations & derivations '''
        
        self.gamma=-0.2*f
        self.sdelta=-1.0/120                             #Stevens & Murphy, 1998
        self.delta= -math.log(2)/0.8
        self.ppr=0.5
        
        
        self.alpha_full= -1.0/30.0
        # This rate specifically applies to full-scale fusion vesicles moving from the endocytosed pool to the RCP (Ryan 1993)  
        
        
        
        if f==1:
            self.lambd= -0.0056
            self.k=-0.00456
            self.beta=-1.0/6.0
            self.alph=-0.00000019
            self.alp=-0.000019
            self.alphpr= -0.004554
            self.x2_low=0.00
            self.x2_high=0.02
            
            self.a=0.65587
            self.b=0.2491
        
        elif f==3:
            self.lambd= -0.0132
            self.k=-0.00571
            self.beta=-1.0/5.0
            self.alph=-0.00000068
            self.alp=-0.000068
            self.alphpr= -0.005697
            
            self.a=0.65
            self.b=0.35
        
        elif f==10:
            self.lambd= -0.0310
            self.k=-0.01489
            self.beta=-1.0/3.0
            self.alph=-0.00000066
            self.alp=-0.000066
            self.alphpr= -0.014811
            self.x2_low=0.05                    #Lower limit on full scale fusion 
            self.x2_high=0.15                   #Upper limit on full scale fusion
            
            self.a=0.3252
            self.b=0.6748
        
        elif f==30:
            self.lambd= -0.0621
            self.k=-0.0538
            self.beta=-1.0/3.0
            self.alph=-0.00000022
            self.alp=-0.000022
            self.alphpr= -0.053356
            self.x2_low=0.25                    #Lower limit on full scale fusion (but why such a value??)
            self.x2_high=0.36                   #Upper limit on full scale fusion
            self.x=10
            
            self.a=0.49920                  #Coefficients of the exponential function used to 
            self.b=0.50070
            
            
        self.str= str(f)+" Hz"             #For the purposes of making documentation easier.
        
        self.x2=ran.choice(GeneralOne.x2range(self))            #Chosen extent of full scale fusion
        #self.x2=0
        
        self.x1= 1 - self.x2                      #Chosen extent of kiss and run
        
        self.epold= ((self.gamma+2*self.sdelta)*self.lambd)/(self.gamma+self.lambd)
        print "Old Epsilon:\t%f" %(self.epold)
        self.ep= (self.lambd/(self.gamma*(self.gamma+ self.lambd)))*((self.gamma)**2 + (self.ppr*self.sdelta)*(self.gamma+self.lambd))
        print "Epsilon:\t%f" %(self.ep)
        self.deprt= -math.log(2)/2.5                        #Chosen rate of departitioning of fm1-43 dye ( based on 2.5 s halftime)
        
        
        
    
obj=GeneralTwo()    
if __name__ == "__main__":
    
    obj.controlpanel()
        