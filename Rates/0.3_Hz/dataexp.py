# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 15:33:06 2018

@author: Koustav
"""

import os
import numpy as np
import math

gamma= 0.06 
# 0.2*0.3
# [ln(2)/(3/(0.2*0.3))], where 3/(0.2*0.3) is the half time (time taken for three vesicles to fuse)
lambd=0.00405   
# EGFP-VAMP Data for 0.3Hz after 100s cutoff: 0.209+0.650*e^(-0.00241x); -0.00241 +/- 0.00015
# EGFP-VAMP Data for 0.3Hz after 50s cutoff: 0.263+0.6424*e^(-0.00305*x); -0.00305 +/- 0.00013
delta=1.0/120 
#Undocking rate obtained from Stevens, Murphy 1998.
def C():
    global gamma
    os.chdir("../../../Data/SVFR Decay Rates")
    svfr=np.genfromtxt("0.3 Hz SVFR Decay.csv", delimiter=",")
    dC(svfr)
    
def dC(svfr):
    global gamma
    C=svfr[:,1]/gamma           #Taking this value to be C
    t=svfr[:,0]
    print t
    dC=[]
    for x in range(0,len(C)-2):
        lc=(C[x+1]-C[x])/(t[x+1]-t[x])
        rc=(C[x+2]-C[x+1])/(t[x+2]-t[x+1])
        mc=(lc+rc)/2
        dC.append(mc)
    dC.extend([0,0])            #Generating dC
    eB=Bgen(t)
    print len(eB)
    #print eB
    dA=[]
    for x in range(0,len(dC)):
        m=dC[x]+svfr[x,1]+delta*C[x]-eB[x]
        dA.append(m)
    print "dA:", dA
    output=np.zeros([len(C),5])
    output[:,0]=t
    output[:,1]=C
    output[:,2]=dC
    output[:,3]=eB
    output[:,4]=dA
    for x in output:
        print output
    os.chdir("../../Code/Rates/0.3_Hz")
    os.mkdir("L_%f_g_%f_alpha_%f_d" %(lambd, gamma, delta))
    os.chdir("L_%f_g_%f_alpha_%f_d" %(lambd, gamma, delta))
    np.savetxt('0.3 Hz Alpha Rate Data.csv', output, delimiter=",", header="0:Time, 1:C, 2:dC/dt, 3:eB, 4:dA/dt",comments="#")
        
def Bgen(t):
        global lambd; global gamma; global delta
        ep= ((gamma+2*delta)*lambd)/(gamma+lambd)
        print "ep:",ep
        b0=24
        eB=[]
        for x in range(0,len(t)):
            b=ep*b0*(0.209+0.650*math.exp(-ep*t[x]))
            eB.append(b)
        return eB
    
C()