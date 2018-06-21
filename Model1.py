# -*- coding: utf-8 -*-
"""
Created on Thu May 24 19:51:37 2018

@author: Koustav
"""

import numpy as np
#import scipy as scp
import random as ran
import math
import matplotlib.pyplot as plt

x2ran=[]
# All Notations for attributes follow as defined from MODEL.
ftrack=30; alphatrack=0
# Keeps track of the numbers of vesicles that haven't fused even once at the last time step.

class Alpha:
    
    def alphameasure(self):
        global alphatrack
        
        for y in range(self.trps):
          if (self.pref[y,3] == 0):
            ch=math.exp(self.phsehalf*self.pref[y,2])
            chec=ran.random()
            if chec > ch:           #Vesicle must be recycled 
               #self.pref[y,2]=0     #Number of times it has recycled gets reset to 0.
               self.tracker[y]=0    #Fluorescence gets reset to 0.
               self.pref[y,3]=1     #Ensure it never meets this filter condition again.
               self.alphak.append(1)
               self.alphak[-1]=len(self.alphak)
               self.alphat.append(self.pref[y,1])
               
    def plotter(self):
        self.alphat.sort()    #WHY AM I MUTHAFUKKING SORTING?
        plt.plot(self.alphat,self.alphak, 'ro')
        plt.xlabel("Time")
        plt.ylabel("Vesicles Leaving RRP")
        plt.savefig("Plot/AlphaKinetics1.png", dpi=200)
        plt.show()
        plt.close()
        print "Madarjhat"
        for x in self.alphak:
            print x,
        print "Time of alpha kinetics:"
        for x in self.alphat:
            print x,

class ThirtyHz(Alpha):
   #EGFP-VAMP data gives decay graph of: 0.1395+0.851*e^(-0.0497*x) || alternatively rate= -0.0642 || 0.1777+0.828*e^(-0.0908*x)
    
    def __init__(self):
        self.xrange=np.array([5,6,7,8,9,10])  #Possible values that RRP size can attain
        self.x=ran.choice(self.xrange)        #Choose one RRP size at random from xrange
        self.trps=30
        self.x2_low=0.38                    #Lower limit on full scale fusion
        self.x2_stepsize=0.01
        self.x2_high=0.51                   #Upper limit on full scale fusion
        self.x2=ran.choice(ThirtyHz.x2range(self))            #Chosen extent of full scale fusion
        self.x1= float(1 - self.x2)                      #Chosen extent of kiss and run
        self.bk= -0.0497                                    #Chosen rate of first fusion
        self.glbk= -0.038                                   #Chosen rate of reuse
        self.deprt= -math.log(2)/2.5                        #Chosen rate of departitioning of fm1-43 dye ( based on 2.5 s halftime)
        self.tracker=np.ones(self.trps)                     #Setting up tracker for individual vesicles with each element containing initial fluorescence.
        #self.trfusion1=np.zeroes(self.trps)                 #Setting up tracker which measures whether an individual vesicle has undergone first fusion
        self.t=0              #Initial time
        self.dt=0.0333          #1/f chosen as time step.
        self.unfsze=0           #Keeps count of number of unfused vesicles
        self.fuskr=0            #Keeps count of number of vesicles that fuse and run.
        self.fulfs=0            #Keeps count of number of vesicles that undergo full fusion.
        self.pref=np.zeros([self.trps,4])  
        '''Keeps track of changes taking place in vesicles, with 1st column storing kind of fusion ("100" for k&r, "200" for full fusion,
        "300" for TRP), 2nd column noting exact time of this occuring and the 3rd one noting the number of such consecutive fusions. '''
        self.pref[:,2]=0
        self.phsehalf= -math.log(2)/8            #Chosen rate of replacement of ageing SV (half-time= 8 cycles)
        self.fsum=[]
        self.ftime=[]
        self.alphak=[]                  #Stores the number of vesicles that have been extinguished at any one point.
        self.alphat=[]
        
    def intpr(self):
        print "\n RRP Size:",self.x," ",
        print "Prob. of Full Scale Fusion:",self.x2," "
        for x in self.tracker:
            print x," ",
        print "\nDetails on all the vesicles at the end:"
        for x in self.pref:
            print x
        
    def firstfusiontracker(self):               #Updates vesicle fluorescence trackers, viz a viz first fluorescence
        global ftrack
        rrpfused= 0.15*self.t*30         # rrpfused documents how many RRP have fused with the membrane FOR THE FIRST TIME.
        retain=0
        if (rrpfused <= self.x):
            retain= (self.trps-self.x)*(0.1395+ 0.851*math.exp((self.bk*self.t))) + self.x- rrpfused
        else:
            retain= (self.trps-self.x)*(0.1395+ 0.851*math.exp((self.bk*self.t)))
        # retain documents how many TRP vesicles haven't undergone first fusion at a particular time.
        print "retain: ", retain
        if( 2>(ftrack-retain) >= 1 ):   #At least one vesicle has fused.
            ftrack=int(math.ceil(retain))
            ind= int(self.trps-ftrack-1) #Keeps track of which element in self.tracker (ie vesicle number) to alter.
            
            if ( ran.random()< self.x2): #Determines what happens if vesicle fusing is full scale fusion
                self.tracker[ind]=0      #All fluorescence discharged.
            else:
                self.tracker[ind] *= math.exp(self.deprt*0.1*ran.choice(range(6,10)))
                #Taking time of kiss and run to be 0.6-0.9, calculating amount of fluorescence discharged and updating accordingly.
        elif(ftrack-retain>=2):
            print "Meep!"  #This is currently a flaw in my code that doesn't make a difference.
                
        else: pass
        
    def counter(self):
        self.fulfs = self.fuskr= self.unfsze=0   #Resetting count to 0.
        for x in self.tracker:
            if x== 1:
                self.unfsze+=1
                print "Meep!"
            elif (1>x>0):
                self.fuskr+=1
            elif x==0:
                self.fulfs+=1
                
    def controlpanel(self):        
        copyt=self.tracker #Preserves a copy of tracker for comparison at T=t
        print copyt
        for x in range(0,6000):
            copyt=tuple(self.tracker)
            self.t+=self.dt
            ThirtyHz.updateparam(self)
            print copyt
            print self.tracker
            for y in range(self.trps):
                if ((0 < self.tracker[y] < copyt[y]) and self.pref[y,2]==0):  #K&R has taken place for  1st time
                    print "truesd"
                    self.pref[y,0]=100                      #Indicates K&R has occurred.
                    self.pref[y,1]=self.t                   #Indicates time instance of fusion.
                    self.pref[y,2]=1                       #Indicates the number of times fusion has taken place.

                elif (copyt[y]>0 and self.tracker[y]==0 and self.pref[y,2]==0):    #Full scale fusion has occured for 1st time
                    print "flaskd"
                    self.pref[y,0]=200                      #Indicates K&R has occurred.
                    self.pref[y,1]=self.t                   #Indicates time instance of fusion.
                    self.pref[y,2]=1                       #Indicates the number of times fusion has taken place.
                   
                else: pass                                  # THIS DOES NOT ACCOMODATE FOR FUSIONS 1st FUSION
                
        ThirtyHz.plotter(self)
        Alpha.plotter(self)
                    
            
    def updateparam(self):
        self.firstfusiontracker()
        self.activatehash()
        self.checkhash()
        
        self.fluores()
        
        
        
    def activatehash(self):             #Go into Narcos mode. Check if time is ripe for vesicle to enter TRP pool (get recycled)
        for y in range(self.trps):
            m=self.pref[y,1]
            if(self.pref[y,0]==100 and (m+0.9-self.t)<self.dt):
                self.pref[y,0]=300
                self.pref[y,1]+=0.9         #Update to reflect time of entering TRP
                self.alphameasure()
            elif(self.pref[y,0]==200 and (m+20-self.t)<self.dt):
                self.pref[y,0]=300
                self.pref[y,1]+=20          #Update to reflect time of entering TRP
                self.alphameasure()
            else: continue
        
    def checkhash(self):                #Check whether a given vesicle remains in TRP or refuses and updates matrices accordingly. 
        for y in range(self.trps):
            if(self.pref[y,0]==300):
               m=self.t - self.pref[y,1]
               ch=math.exp((self.glbk*m))
               if( ran.random() > ch):              #Re-fusion has taken place.
                   self.pref[y,1]= self.t
                   self.pref[y,2]+=1
                   print "Ganfu!"
                   chec=ran.random()
                   if (chec <= self.x2):           #Full scale re-fusion has taken place.
                       self.pref[y,0]=200
                       self.tracker[y]=0                #Setting fluorescence parameter to be 0
                   elif (chec > self.x2):            # K&R re-fusion has taken place.
                       self.pref[y,0]=100
                       self.tracker[y]*= math.exp(self.deprt*0.1*ran.choice(range(6,10)))
                       '''Taking time of kiss and run to be 0.6-0.9, calculating amount of fluorescence
                       discharged and updating accordingly.'''
               else: continue
           

  
    def fluores(self):
        f=0
        for y in range(self.trps):
            f+=self.tracker[y]
        self.fsum.append(f)
        self.ftime.append(self.t)
        
    def plotter(self):
        plt.plot(self.ftime,self.fsum, 'c-')
        plt.xlabel("Time")
        plt.ylabel("Flourescence")
        plt.savefig("Plot/Destaining1.png", dpi=200)
        plt.show()
        plt.close()
    
    def x2range(self):
        x= self.x2_low
        global x2ran
        while ( x<=self.x2_high ):
            x2ran.append(x)
            x+=self.x2_stepsize
        else:
            return x2ran
        
        
obj=ThirtyHz()
if __name__=="__main__":
    obj.intpr()
    obj.controlpanel()        
    obj.intpr()
    