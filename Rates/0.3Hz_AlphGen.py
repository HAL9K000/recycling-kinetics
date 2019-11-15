# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 14:09:33 2019

@author: Koustav
"""

import math

def alpha_gen():
    
    k= 0.00408
    lambd=0.00406
    sdelt=1/120.0
    gamma=0.06
    ep= lambd*((gamma**2)+sdelt*(gamma+lambd))/(gamma*(gamma+lambd))
    print "Epsilon:\t %6.5f" %(ep)
    delt=math.log(2)/0.9 #Assuming time constant of 0.9s
    delt_pr= math.log(2)*((1/0.9)-(1/2.5)) #Chosen rate of departitioning of fm1-43 dye ( based on 2.5 s halftime)
    
    coef2= ((lambd*sdelt)/(gamma*(k+ep)+lambd))-1
    
    alp= k*(1+ (delt_pr*coef2/(k-delt)))
    
    print "Alpha Prime Value :\t %6.5f" %(alp)
    
    
alpha_gen()
    