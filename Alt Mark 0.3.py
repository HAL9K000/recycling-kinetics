# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 12:56:36 2018

@author: Koustav
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import os
import random as ran
from scipy.optimize import curve_fit

eptrack=0

class PointthreeHz:
# From EGFP-VAMP Data, we have, best fit: 0.3191 + 0.6456e^(-0.0406*x)   
    def __init__(self):
        global eptrack
        self.xrange=np.array([5,6,7,8,9,10])  #Possible values that RRP size can attain
        self.x=ran.choice(self.xrange)        #Choose one RRP size at random from xrange
        self.x=10
        self.trps=30                        #Size of TRP
        self.x2_low=0.10                    #Lower limit on full scale fusion (but why such a value??)
        self.x2_stepsize=0.01
        self.x2_high=0.19                   #Upper limit on full scale fusion
        self.x2=ran.choice(PointthreeHz.x2range(self))            #Chosen extent of full scale fusion
        self.x1= 1 - self.x2                      #Chosen extent of kiss and run
        self.gamma= -0.06                                    #Look up the model for meanings of specific parameters
        self.delta=-1.0/120                             #Stevens & Murphy, 2018
        self.lambd= -0.00406
        self.alphpr= -0.0099
        self.ep= ((self.gamma+2*self.delta)*self.lambd)/(self.gamma+self.lambd)
        self.deprt= -math.log(2)/2.5                        #Chosen rate of departitioning of fm1-43 dye ( based on 2.5 s halftime)
        self.tracker=np.ones(self.trps)                     #Setting up tracker for individual vesicles with each element containing initial fluorescence.
        self.t=0              #Initial time
        self.dt=1          # Chosen as time step.
        self.pref=np.zeros([self.trps,4])
        '''Keeps track of changes taking place in vesicles, with 1st column storing kind of fusion 
        ("100" for k&r, "200" for full fusion, "300" for endocytosed, "400" for RRP, "0" for RCP unfused), 2nd column 
        noting exact time of this occuring and the 3rd one noting the number of such consecutive fusions.'''
        self.pref[0:self.x,0]=400
        eptrack=self.trps-self.x                    #Notes number of vesicles in RCP at the beginning
        self.fsum=[]                        #Measures total fluorescence at each time step.
        self.rrpsiz=[]
        self.endsiz=[]
        self.ftime=[]
        
    def controlpanel(self):
        print "x :=", self.x, "\tx2 :=", self.x2 ,"\nInitially:"
        for x in range(self.trps):
            print self.pref[x]
        for x in range(0,800):
            self.updateparam()
            self.t+=self.dt
        for x in self.pref:
            print x
        PointthreeHz.accountancy(self)
        
    def updateparam(self):
        self.rcp_rrp()
        self.fusion()
        self.endo()
        self.priming()
        self.fluores()
        
    def rcp_rrp(self):                     #Determines how many vesicles have moved out from RCP to RRP based on kinetic data.
        global eptrack
        retain= (self.trps-self.x)*(0.88+ 0.12*math.exp(self.ep*self.t))
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
                
                ch=0.6895+ 0.3105*math.exp(self.ep*m)
                if (ran.random()>ch):
                    #Vesicle is primed and docked at RRP
                    self.pref[y,0]=400
                    self.pref[y,1]=self.t   #Updating time of entry to RRP.
        
    def fusion(self):                   #Goes into Narcos mode, checking for and updating vesicle parameters as they fuse with membrane.
        for y in range(self.trps):
            if (self.pref[y,0]==400):       #Checkes whether a vesicle is in the RRP.
                m=self.t - self.pref[y,1]        #Notes the amount of time a vesicle has spent in the RRP
                ch=math.exp((self.delta*m))
                if( ran.random() > ch):
                    #Vesicle is undocked from RRP to RCP
                    self.pref[y,0]=0
                    self.pref[y,1]=self.t   #Updating time of entry to RRP.
                    print "BC"
                    continue #Move to next iteration
                    print "CD"
                ch=math.exp((self.gamma*m))
                if( ran.random() > ch):              #Fusion has taken place.
                   print "Ganfu!"
                   if (self.pref[y,2]==0):              #Fusion for the first time.
                        print "Meep"
                        chec=ran.random()
                        if (chec <= self.x2):           #Full scale fusion has taken place.
                            print "Bo"
                            self.pref[y,0]=200
                            self.tracker[y]=0                #Setting fluorescence parameter to be 0
                        elif (chec > self.x2):            # K&R fusion has taken place.
                            self.pref[y,0]=100
                            self.tracker[y]*= math.exp(self.deprt*0.1*ran.choice(range(6,10)))
                            '''Taking time of kiss and run to be 0.6-0.9, calculating amount of fluorescence
                            discharged and updating accordingly.'''
                   else:                                #Refusion has taken place
                        self.pref[y,0]=100
                        self.tracker[y]*= math.exp(self.deprt*0.1*ran.choice(range(6,10)))
                        
                   self.pref[y,1]= self.t     #Update the time parameter to reflect fusion
                   self.pref[y,2]+=1        #Notes the number of times fusion has taken place.
            else: continue
        
    def endo(self):  #Checks whether vesicle is endocytosed and updates parameters accordingly.
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
        
    def priming(self):  #Checks whether endocytosed vesicle is primed and docked and updates parameters accordingly.
        for y in range(self.trps):
            if (self.pref[y,0]==300):
                #Vesicle is in endocytosed compartment.
                m=self.t-self.pref[y,1]         #Notes the time spent by vesicle as endocytosed.
                ch=math.exp((self.alphpr*m))
                if (ran.random()>ch):
                    #Vesicle is primed and docked at RRp
                    self.pref[y,0]=400
                    self.pref[y,1]=self.t   #Updating time of entry to RRP.
            
    def fluores(self):
        f=0; rrp=0 ; end=0
        for y in range(self.trps):
            f+=self.tracker[y]
            if (self.pref[y,0]==400):               #Measures RRP size at each step
                rrp+=1
            elif (self.pref[y,0]==300):
                end+=1
        self.rrpsiz.append(rrp)
        self.endsiz.append(end)
        self.fsum.append(f)
        self.ftime.append(self.t)
        
    def f(self, x,a,b,c):             #Used by curve-fitting algorithm
        return a+ b*np.exp(-c*x)
        
    def accountancy(self):          #Plots, analyses and saves various generated results.
        os.chdir("Plot/0.3_Hz")
        fm143=np.genfromtxt("0.3 Hz FM143.csv", delimiter=',')
        fm143[:,1]*=30
        plt.plot(fm143[:,0] ,fm143[:,1] , 'bo', markerfacecolor='none', markeredgecolor='b')
        pop, pco = curve_fit(self.f, fm143[:,0], fm143[:,1], bounds=([15.0,0,0],[30.0,15.0, 0.05]))
        plt.plot(fm143[:,0], self.f( fm143[:,0], *pop), 'c--', label='FM Fit: %5.4f + %5.4f*e^(-%5.4fx)\n' % tuple(pop))
        #Plotting FM1-43 data.
        a=-self.alphpr; g=-self.gamma; d=-self.delta
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
        f.write("Alpha: %f,\t Gamma: %f,\t Lambda: %f\t Delta: %f\n" %(self.alphpr, self.gamma, self.lambd, -d))
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
    
    def x2range(self):
        x=[]; i=self.x2_low
        while(i<= self.x2_high):
            x.append(i)
            i+=self.x2_stepsize
        return x
        
    
obj=PointthreeHz()
if __name__ == '__main__':
    obj.controlpanel()
