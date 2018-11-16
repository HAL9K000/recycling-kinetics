# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 13:39:22 2018

@author: Koustav
"""

import math
import os
'''
Goal is to use alpha prime (Apr) values from 0.3 Hz model and fit them into the polynomial equation to determine alpha for higher 
frequency stimulations
'''
alpr=beta=gamma=delta=k=0.0
print "Enter Frequency Number: \t",
f=float(raw_input())
if( os.path.isdir("%2.1f_Hz" %(f))==False):
    os.mkdir("%2.1f_Hz" %(f))
os.chdir("%2.1f_Hz" %(f))

def init():
    global alpr; global beta; global gamma; global delta; global k; global f
    alpr=float(raw_input("Enter Alpha Prime Value: \t"))
    alpr= 0.0025
    #beta=float(raw_input("Enter Beta Value: \t"))
    beta=1.0/5.0
    gamma=0.2*f
    #delta= float(raw_input("Enter Undocking Rate Value: \t"))
    delta=1.0/120.0
    k= float(raw_input("Enter K Value: \t"))
    rootfinder()

def rootfinder():
    global alpr; global beta; global gamma; global delta; global k; global f
    
    b= -(beta*(1+(gamma/k))+gamma+delta- (gamma*alpr/k))
    c= gamma*beta
    r1=(-b+math.sqrt(b*b-4*c))/(2.0)
    r2=-(b+math.sqrt(b*b-4*c))/(2.0)
    a1=r1-alpr; a2=r2-alpr
    
    if (os.path.isdir("B_%f_g_%f_a'_%f_k_%f_d_%f_" %(beta, gamma, alpr, k, delta))== False):
        os.mkdir("B_%f_g_%f_a'_%f_k_%f_d_%f_" %(beta, gamma, alpr, k, delta))
    os.chdir("B_%f_g_%f_a'_%f_k_%f_d_%f_" %(beta, gamma, alpr, k, delta))
    
    f=open("root_data.txt", 'w')
    f.write("Possible Roots:\n Alpha=%5.4f;\t Alpha=%5.4f\n" %(a1,a2))
    f.write("Other Data: \n Alpr: %f,\t Beta=%f,\t Gamma: %f,\t Delta: %f,\t K: %f\n" %(alpr, beta, gamma, delta, k))
    f.flush(); f.close()
    print( "Possible Roots:\n Alpha=%5.4f;\t Alpha=%5.4f\n" %(a1,a2))

init()
    
