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
'''
SVFR Data (Best Fit):
    0.3 Hz: 0.687+0.125*e^(-0.00768*x)
    1   Hz: 0.476+0.239*e^(-0.00456*x)
    3   Hz: 0.420+0.302*e^(-0.00571*x)
    10  Hz: 0.178+0.370*e^(-0.01489*x)
    30  Hz: 0.237+0.246*e^(-0.0538*x)
    General Relation b/w half-time(T) and f: T= 9.9073+ 166.718*e^(-0.1452*f)
'''
class PointthreeHz:
# From EGFP-VAMP Data, we have, best fit: 0.3191 + 0.6456e^(-0.0406*x)   
    def __init__(self):
        #global eptrack
        self.xrange=np.array([5,6,7,8,9,10])  #Possible values that RRP size can attain
        self.x=ran.choice(self.xrange)        #Choose one RRP size at random from xrange
        
        self.ratepicker()
        
        self.tracker=np.ones(self.trps)                     #Setting up tracker for individual vesicles with each element containing initial fluorescence.
        self.t=0              #Initial time
        self.dt=1          # Chosen as time step.
        self.pref=np.zeros([self.trps,4])
        '''Keeps track of changes taking place in vesicles, with 1st column storing kind of fusion 
        ("100" for K&R, "200" for Full Fusion, "300" for Endocytosed ("301" for K&R, "302" for full fusion), "400" for RRP, 
        "0" for RCP unfused, "500" for ARCP, 2nd column noting exact time of this occuring,
        the 3rd one noting the number of such consecutive fusions,while the 4th stores a random number
         that determines the probability of transfer of each vesicle from one to another pool'''
        self.pref[0:self.x,0]=400       #
        #eptrack=self.trps-self.x                    #Notes number of vesicles in RCP at the beginning
        self.fsum=[]                        #Measures total fluorescence at each time step.
        self.rrpsiz=[]
        self.endsiz=[]
        self.ftime=[]
        
        for r in range(0,self.x):
            self.pref[r,3]= ran.random() #Generating random number for these beauts.
            
        for r in range(self.x,self.trps):
            self.pref[r,3]= ran.random() #Generating random number for these beauts
        
    def controlpanel(self):
        print "x :=", self.x, "\tx2 :=", self.x2 ,"\nInitially:"
        for x in range(self.trps):
            print self.pref[x]
        for x in range(0,800):
            self.updateparam()
            self.t+=self.dt
        for x in self.pref:
            print x
        self.accountancy()
        
    def updateparam(self):
        self.rcp_rrp()
        self.fusion()
        self.endo()
        self.priming()
        self.fluores()
        
    def rcp_rrp(self):                     #Determines how many vesicles have moved out from RCP to RRP based on kinetic data.
        '''global eptrack
        retain= (self.trps-self.x)*(0.93+ 0.07*math.exp(self.ep*self.t))
        #print "Retain:", retain
        # retain documents how many TRP vesicles haven't made the rcp-rrp transit at a particular time.
        if (eptrack-retain >=1):
            print "Kookoo"
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
                
                ch=0.93+ 0.07*math.exp(self.ep*m)
                if (ran.random()>ch):
                    #Vesicle is primed and docked at RRP
                    self.pref[y,0]=400
                    self.pref[y,1]=self.t   #Updating time of entry to RRP.'''
                    
                    
        #Old Code above.
        
        for y in range(self.trps):
            
            if (self.pref[y,0]==0 and self.pref[y,1] == 0):
                #Vesicle has never left Latent RCP even once.
                #print "Saheb Bibi Golam"
                
                m=self.t #Time spent by Vesicle in Latent RCP.
                
                ch=self.a+ self.b*math.exp(self.ep*self.t)
                if (self.pref[y,3]>ch):
                    #Vesicle is primed and docked at RRP
                    self.pref[y,0]=400
                    self.pref[y,1]=self.t   #Updating time of entry to RRP.
                    self.pref[y,3]=ran.random()
            
            if (self.pref[y,0]==0 and self.pref[y,1]>0):
                #Vesicle has been undocked to  Latent RCP.
                m=self.t-self.pref[y,1]         #Notes the time spent by vesicle in RCP.
                
                ch=self.a+ self.b*math.exp(self.ep*m)
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
                if( self.pref[y,3] > ch):
                    #Vesicle is undocked from RRP to LRCP
                    self.pref[y,0]=0
                    self.pref[y,1]=self.t   #Updating time of entry to RCP.
                    self.pref[y,3]= ran.random()
                    print "BC"
                    continue #Move to next iteration
                    
                ch=math.exp((self.gamma*m))
                if( self.pref[y,3] > ch):              #Fusion has taken place.
                   print "LoFi"
                   if (self.pref[y,2]==0):              #Fusion for the first time.
                        print "Meep"
                        chec=ran.random()
                        if (chec <= self.x2):           #Full scale fusion has taken place.
                            print "Impossible"
                            self.pref[y,0]=200
                            self.tracker[y]=0                #Setting fluorescence parameter to be 0
                            self.pref[y,3]=ran.uniform(0.3, 0.6)      
                            #Updating random-number stored so that the vesicle will undergo endocytosis between 15-35 s
                        elif (chec > self.x2):            # K&R fusion has taken place.
                            self.pref[y,0]=100
                            #self.tracker[y]*= math.exp(self.deprt*0.1*9)
                            self.pref[y,3]=ran.uniform(0.32, 0.68) 
                            #K&R vesicle  must not spend more than (1.31) 2 seconds docked at the membrane or less than 0.45 s
                            # With a mean of 0.8 seconds docked at the membrane
                            t_calc= (math.log(self.pref[y,3])/self.delta)
                            #t_calc is the amount of time spent docked at the membrane.
                        
                            self.tracker[y]*= math.exp(self.deprt*t_calc)
                   else:                                #Refusion has taken place (Exclusively K&R)
                        self.pref[y,0]=100
                        self.pref[y,3]=ran.uniform(0.32, 0.68) 
                        #K&R vesicle  must not spend more than (1.31) 2 seconds docked at the membrane or less than 0.45 s
                        # With a mean of 0.8 seconds docked at the membrane
                        t_calc= (math.log(self.pref[y,3])/self.delta)
                        
                        self.tracker[y]*= math.exp(self.deprt*t_calc)
                        
                   self.pref[y,1]= self.t     #Update the time parameter to reflect fusion
                   self.pref[y,2]+=1        #Notes the number of times fusion has taken place.
            else: continue
        
    def endo(self):  #Checks whether vesicle is endocytosed and updates parameters accordingly.
        for y in range(self.trps):
            m=self.t - self.pref[y,1]
            ch1=math.exp((self.delta*m))
            ch2=math.exp(((math.log(2)/20.0)*m))
            if( self.pref[y,0]==100 and self.pref[y,3] > ch1):
                #K&R vesicle is endocytosed.
                self.pref[y,0]=300          
                self.pref[y,1]+=m         #Time of entering endocytosed compartment is updated.
                self.pref[y,3]= ran.random()
            elif (self.pref[y,0]==200 and self.pref[y,3] > ch2):
                #Full fusion vesicles are endocytosed
                self.pref[y,0]=300          
                self.pref[y,1]+=m         #Time of entering endocytosed compartment is updated.
                self.pref[y,3]= ran.random()
        
    def priming(self):  #Checks whether endocytosed vesicle is primed and docked and updates parameters accordingly.
        for y in range(self.trps):
            if (self.pref[y,0]==300):
                #Vesicle is in endocytosed compartment.
                m=self.t-self.pref[y,1]         #Notes the time spent by vesicle as endocytosed.
                ch=math.exp((self.alphpr*m))
                if (self.pref[y,3]>ch):
                    #Vesicle is primed and docked at RRp
                    self.pref[y,0]=400
                    self.pref[y,1]=self.t   #Updating time of entry to RRP.
                    self.pref[y,3]= ran.random()
            
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
        os.chdir("Plot/0.3 Hz")
        fm143=np.genfromtxt("0.3 Hz FM143.csv", delimiter=',')
        fm143[:,1]*=30
        plt.plot(fm143[:,0] ,fm143[:,1] , 'bo', markerfacecolor='none', markeredgecolor='b')
        pop, pco = curve_fit(self.f, fm143[:,0], fm143[:,1], bounds=([15.0,0,0],[30.0,15.0, 0.05]))
        plt.plot(fm143[:,0], self.f( fm143[:,0], *pop), 'c--', label='FM Fit: %5.4f + %5.4f*e^(-%5.4fx)\n' % tuple(pop))
        #Plotting FM1-43 data.
        a=-self.alphpr; g=-self.gamma; d=-self.sdelta
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
        
        
    def ratepicker(self):
        #Set all parameters appropriately.
        
        self.x2_low=0.00                    #Lower limit on full scale fusion (but why such a value??)
        self.x2_stepsize=0.01
        self.x2_high=0.00                  #Upper limit on full scale fusion
        
        self.x=6
        self.trps=30                        #Size of TRP
        
        self.x2=ran.choice(PointthreeHz.x2range(self))            #Chosen extent of full scale fusion
        self.x2= 0.00       #NO FULL FUSION.
        self.x1= 1 - self.x2                      #ONLY K&R.
        
        
        self.gamma= -0.06                                    #Look up the model for meanings of specific parameters
        self.sdelta=-1.0/120                             #Stevens & Murphy, 2018
        self.delta=-math.log(2)/0.8                  #Assuming time constant of 0.8s to perform K&R.
        self.delta_pr= -math.log(2)*((1/0.8)-(1/2.5))            #Amount of 
        self.lambd= -0.00406
        self.alphpr= -0.00684
        
        self.a=0.875
        self.b=0.125     #Coefficients of release from LRCP.
        
        self.ppr=1.0
        
        
        self.epold= ((self.gamma+2*self.sdelta)*self.lambd)/(self.gamma+self.lambd)
        print "Old Epsilon:\t%f" %(self.epold)
        self.ep= (self.lambd/(self.gamma*(self.gamma+ self.lambd)))*((self.gamma)**2 + (self.ppr*self.sdelta)*(self.gamma+self.lambd))
        self.ep= self.lambd*((self.gamma**2)+self.sdelta*(self.gamma+self.lambd))/(self.gamma*(self.gamma+self.lambd))
        print "Epsilon:\t%f" %(self.ep)
        self.deprt= -math.log(2)/2.5                        #Chosen rate of departitioning of fm1-43 dye ( based on 2.5 s halftime)
        
        
        
        
    
    def x2range(self):
        x=[]; i=self.x2_low
        while(i<= self.x2_high):
            x.append(i)
            i+=self.x2_stepsize
        return x
        
    
obj=PointthreeHz()
if __name__ == '__main__':
    obj.controlpanel()
