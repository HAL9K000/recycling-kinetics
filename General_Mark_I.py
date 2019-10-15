# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 12:47:53 2018

@author: Koustav
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import os
import random as ran
from scipy.optimize import curve_fit

eptrack=0


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

class GeneralOne:
    
    def __init__(self):
        global eptrack
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
        eptrack=self.trps-self.x                    #Notes number of vesicles in RCP at the beginning
        self.fsum=[]                        #Measures total fluorescence at each time step.
        self.rrpsiz=[]
        self.endsiz=[]
        self.ftime=[]
        self.str="10 Hz"             #For the purposes of making documentation easier.
        
    def controlpanel(self):
        print "RRP Size :=", self.x, "\tx2 :=", self.x2 ,"\nInitially.\n"
        for x in range(self.trps):
            print self.pref[x]
        for x in range(0,800):
            self.updateparam()
            self.t+=self.dt
        for x in self.pref:
            print x
        GeneralOne.accountancy(self)
        
    def updateparam(self):
        self.lrcp_rrp()
        self.fusion()
        self.endo()
        self.priming()
        self.acrp_rrp()
        self.fluores()
        
    def lrcp_rrp(self):                     #Determines how many vesicles have moved out from unfused RCP to RRP based on kinetic data.
        global eptrack
        retain= (self.trps-self.x)*(0.3252+ 0.6748*math.exp(self.ep*self.t))
        #print "Retain:", retain
        # retain documents how many TRP vesicles haven't made the rcp-rrp transit at a particular time.
        if (eptrack-retain >=1):
            m=self.trps-eptrack         #Notes the vesicle number that will undergo transit
            k=int(eptrack-math.ceil(retain))        #In the event that there are more than one vesicles making the transit.
            print "k:",k
            for i in range(0,k):
                a=int(m+i)
                print "a:",a
                self.pref[a,0]=400          #Change tag of vesicle to RRP
                self.pref[a,1]=self.t       #Change timing to reflect entry to RRP
            eptrack=math.ceil(retain)       #Update eptrack to reflect current number of unfused vesicles
            for x in range(0,self.trps):
                print self.pref[x,0],
                
        for y in range(self.trps):
            if (self.pref[y,0]==0 and self.pref[y,1]>0):
                #Vesicle has been undocked to RCP.
                m=self.t-self.pref[y,1]         #Notes the time spent by vesicle in RCP.
                
                ch=0.3252+ 0.6748*math.exp(self.ep*self.t)
                if (ran.random()>ch):
                    #Vesicle is primed and docked at RRP
                    self.pref[y,0]=400
                    self.pref[y,1]=self.t   #Updating time of entry to RRP.
        
    def fusion(self):                   #Goes into Narcos mode, checking for and updating vesicle parameters as they fuse with membrane.
        for y in range(self.trps):
            if (self.pref[y,0]==400):       #Checkes whether a vesicle is in the RRP.
                m=self.t - self.pref[y,1]        #Notes the amount of time a vesicle has spent in the RRP
                ch=math.exp((self.delta*m))
                #ch1=math.exp((self.deltap*m))
                if( ran.random() > ch):
                    #Vesicle is undocked from RRP......
                    
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
                '''if(ran.random()>ch1):
                    #Vesicle is undocked from RRP to Active RCP
                    self.pref[y,0]=500
                    self.pref[y,1]=self.t   #Updating time of entry to Active RCP.
                    print "CB"
                    continue #Move to next iteration'''
                    
                ch=math.exp((self.gamma*m))
                if( ran.random() > ch):              #Fusion has taken place.
                   print "CD!"
                   chec=ran.random()
                   if (chec <= self.x2):           #Full scale fusion has taken place.
                       print "Bo"
                       self.pref[y,0]=200
                       self.tracker[y]=0                #Setting fluorescence parameter to be 0
                   elif (chec > self.x2):            # K&R fusion has taken place.
                        self.pref[y,0]=100
                        #self.tracker[y]*= math.exp(self.deprt*0.1*9)
                        self.tracker[y]*= math.exp(self.deprt*0.1*ran.choice(range(6,10)))
                        '''Taking time of kiss and run to be 0.6-0.9, calculating amount of fluorescence
                        discharged and updating accordingly.'''
                        
                   self.pref[y,1]= self.t     #Update the time parameter to reflect fusion
                   self.pref[y,2]+=1        #Notes the number of times fusion has taken place.
            else: continue
    
    def endo(self):   #Checks whether vesicle is endocytosed and updates parameters accordingly.
        for y in range(self.trps):
            m=self.pref[y,1]
            if( self.pref[y,0]==100 and (m+0.9-self.t)<self.dt):
                #K&R vesicle is endocytosed.
                self.pref[y,0]=300          
                self.pref[y,1]+=0.9         #Time of entering endocytosed compartment is updated.
            elif (self.pref[y,0]==200 and (m+20-self.t)<self.dt):
                #Full fusion vesicles are endocytosed
                self.pref[y,0]=300          
                self.pref[y,1]+=20         #Time of entering endocytosed compartment is updated.
                
    def priming(self):  #Checks whether endocytosed vesicle is primed and updates parameters accordingly.
        for y in range(self.trps):
            if (self.pref[y,0]==300):
                #Vesicle is in endocytosed compartment.
                m=self.t-self.pref[y,1]         #Notes the time spent by vesicle as endocytosed.
                ch=math.exp((self.alphpr*m))
                if (ran.random()>ch):
                    #Vesicle is primed and docked at RRp
                    self.pref[y,0]=400
                    self.pref[y,1]=self.t   #Updating time of entry to RRP.
                    continue
                ch=math.exp((self.alph*m))
                if (ran.random()>ch):
                    #Vesicle is primed and docked at ARCP
                    self.pref[y,0]=500
                    self.pref[y,1]=self.t   #Updating time of entry to ARCP.
                    
    def acrp_rrp(self):  #Checks whether vesicle in Act RCP gets docked at the RRP
        for y in range(self.trps):
            if (self.pref[y,0]==500):
                #Vesicle is in Active RCP
                m=self.t- self.pref[y,1]
                ch=math.exp(self.beta*m)
                if (ran.random()> ch):
                    #Vesicle gets docked at RRP
                    self.pref[y,0]=400
                    self.pref[y,1]=self.t #Updating time of entry to RRP.
           
    
    def fluores(self):
        f=0; rrp=0 ; end=0
        for y in range(self.trps):
            f+=self.tracker[y]
            if (self.pref[y,0]==400):               #Measures RRP size at each step
                rrp+=1
            elif (self.pref[y,0]==300):             #Measures endocytosed pool size at each step.
                end+=1
        self.rrpsiz.append(rrp)
        self.endsiz.append(end)
        self.fsum.append(f)
        self.ftime.append(self.t)
        
    def f(self, x,a,b,c):             #Used by curve-fitting algorithm
        return a+ b*np.exp(-c*x)
    
    def accountancy(self):          #Plots, analyses and saves various generated results.
        os.chdir("Plot")
        if( os.path.isdir(self.str)==False):
            os.mkdir(self.str)
        os.chdir(self.str)
        fm143=np.genfromtxt("%s FM143.csv" %(self.str), delimiter=',')
        fm143[:,1]*=30
        plt.plot(fm143[:,0] ,fm143[:,1] , 'bo', markerfacecolor='none', markeredgecolor='b')
        pop, pco = curve_fit(self.f, fm143[:,0], fm143[:,1], bounds=([10.0,0,0],[30.0,20.0, 0.05]))
        plt.plot(fm143[:,0], self.f( fm143[:,0], *pop), 'c--', label='FM Fit: %5.4f + %5.4f*e^(-%5.4fx)\n' % tuple(pop))
        #Plotting FM1-43 data.
        a=-self.alph; g=-self.gamma; d=-self.delta
        if (os.path.isdir("t_%f_t" %(self.dt))==False):
            os.mkdir("t_%f_t" %(self.dt))
        #Making various directories to store results, if they do not exist to begin with.
        os.chdir("t_%f_t" %(self.dt))
        if (os.path.isdir("alp_%f_g_%f_d_%f" %(a, g, d))==False):
            os.mkdir("alp_%f_g_%f_d_%f" %(a, g, d)) 
        os.chdir("alp_%f_g_%f_d_%f" %(a, g, d))
        
        if (os.path.isdir("rrp_%d_prob_%f" %(self.x, self.x2))==False):
            os.mkdir("rrp_%d_prob_%f" %(self.x, self.x2)) 
        os.chdir("rrp_%d_prob_%f" %(self.x, self.x2))
        
        
        plt.plot(self.ftime , self.fsum, 'ro', markerfacecolor='none', markeredgecolor='r')
        popt, pcov = curve_fit(self.f, self.ftime, self.fsum)
        xdata=np.array(self.ftime)
        plt.plot(xdata, self.f(xdata, *popt), 'm--', label='Th Fit: %5.4f + %5.4f*e^(-%5.4fx)' % tuple(popt))
        plt.xlabel("Time")
        plt.ylim(0,30)
        plt.ylabel("Fluorescence")
        plt.legend()
        plt.savefig("Destaining.png", dpi=200)
        plt.show()
        plt.close()
        
        self.addtnplot()
        # Perform additional plots
        f=open("log.txt", 'w')
        f.write("RRP Size: %d,\t Full Scale Fusion Scale: %f\n" %(self.x, self.x2))
        f.write("Alpha: %f,\t Gamma: %f,\t Lambda: %f\t Delta: %f\n" %(self.alph, self.gamma, self.lambd, -d))
        f.write("Theoretical Fit: %5.4f + %5.4f*e^(-%5.5fx)\n" % tuple(popt))
        f.write("Optimum Fit: %5.4f + %5.4f*e^(-%5.5fx)" % tuple(pop))
        f.flush(); f.close() #Writing home some key parameter data to a file
        
    def addtnplot(self):
        plt.plot(self.ftime , self.rrpsiz, 'ko', markerfacecolor='none', markeredgecolor='k')
        plt.ylabel("RRP Size")
        plt.xlabel("Time")
        plt.savefig("RRP.png",dpi=200)
        plt.show()
        plt.close()
        plt.plot(self.ftime , self.endsiz, 'mo', markerfacecolor='none', markeredgecolor='m')
        plt.ylabel("Endocytosed Pool Size")
        plt.xlabel("Time")
        plt.savefig("Endo.png",dpi=200)
        plt.show()
        plt.close()
        
        
    def x2range(self):                      #Calculates extent of full scale fusion b/w two hard limits.
        x=[]; i=self.x2_low
        while(i<= self.x2_high):
            x.append(i)
            i+=self.x2_stepsize
        return x
    

obj=GeneralOne()    
if __name__ == "__main__":
    
    obj.controlpanel()