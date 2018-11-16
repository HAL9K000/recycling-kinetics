# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 16:49:28 2018

@author: Koustav
"""

import numpy as np
import os
import math
#  SVFR decay fit for 0.3 Hz: 0.687 + 0.125e^(-0.00769)
#  SVFR decay fit for 1.0 Hz: 0.476 + 0.239*e^(-0.00456)
# Predicted half-time from SVFR exponential decay: 172.52 ( r= -0.0040177)

def manufacture_svfr():
    m=0
    a= 0.687; b= 0.125; c= -0.00401
    t=[]
    svfr=[]
    for x in range(0,35):
        s=x*600.0/35.0
        t.append(s)
    for x in t:        
        svfr.append([x, a+b*math.exp(c*x)])
    os.chdir("../../Data/SVFR Decay Rates")
    np.savetxt('Predicted 0.3 Hz SVFR Decay.csv', svfr, delimiter=",", header="0:Time, 1:SVFR",comments="#")
    
manufacture_svfr()