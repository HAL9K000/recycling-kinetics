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
alpr=beta=gamma=sdelta=k=lambd=0.0

qpr=0.0 #Percentage full-fusion
eta=0.0 #Rate of fluoroscence loss
deltapr=0.0
ep=0.0

'''print "Enter Frequency Number: \t",
f=float(raw_input())
if( os.path.isdir("%2.1f_Hz" %(f))==False):
    os.mkdir("%2.1f_Hz" %(f))
os.chdir("%2.1f_Hz" %(f))'''

def init():
    global alpr; global beta; global gamma; global sdelta; global k; global f; global p; global qpr; global eta; global deltapr
    global lambd; global ep
    #alpr= 0.0025
    #beta=float(raw_input("Enter Beta Value: \t"))
    beta=-1.0/5.0
    sdelta=-1.0/120.0
    p=0.5
    
    f= float(raw_input('Enter input frequency:\t'))
    
    
    k,qpr,beta, lambd= assigner(f, k, qpr, beta, lambd)
    
    
    eta= -((math.log(2)/0.9)*(1-qpr) +(math.log(2)/2.5)*qpr)
    
    deltapr=-math.log(2)*((1/0.9)-(1/2.5))*(1-qpr)
    
    
    gamma= -0.2*f
    #delta= float(raw_input("Enter Undocking Rate Value: \t"))
    
    print "K:\t%f" %(k)
    print "Lambda:\t%f" %(lambd)
    print "Gamma:\t%f" %(gamma)
    
    #k= float(raw_input("Enter K Value: \t"))
    rootfinder3()
    
def assigner(f,k,qpr,beta, lambd):
    
    if f==1:
        lambd= -0.0056
        k=-0.00456
        beta=-1.0/6.0
        qpr=0.05
    elif f==3:
        lambd= -0.0132
        k=-0.00571
        beta=-1.0/5.0
        qpr=0.13
    elif f==10:
        lambd= -0.0310
        k=-0.01489
        beta=-1.0/3.0
        qpr=0.13
    elif f==30:
        lambd= -0.0621
        k=-0.0538
        beta=-1.0/3.0
        qpr=0.3
        
    return k,qpr,beta, lambd

def rootfinder1():
    global alpr; global beta; global gamma; global sdelta; 
    global deltapr; global eta; global p; global qpr;
    global k; global f; global ep
    
    
    q=1-qpr; ppr=1-p;
    
    ep= lambd*((gamma**2)+(sdelta*gamma*ppr)+(lambd*ppr*sdelta))/(gamma*(gamma+lambd))
    print "Epsilon Value:\t%f" %(ep)
    
    theta= gamma +sdelta*(1+(lambd*p)/gamma)
    b_dot= theta - ep*(1+gamma/lambd)
    
    print "Theta:\t%f" %(theta)
    print "B_dot:\t%f" %(b_dot)
    
    r= ppr*sdelta*beta*b_dot/(b_dot*(beta-k)-ppr*beta*sdelta)
    
    print "R:\t%f" %(r)
    
    num= deltapr*q*(((r*b_dot)/(1+r))+k-theta) +eta-k
    denom= (deltapr*q*gamma/k) -(1 - eta/k)*(b_dot/(1+r))
    
    alpr= num/denom
    
    print "Alpha Prime Value:\t%f" %(alpr)
    
    
def rootfinder2(): #Fixed errors, strict segregation of K&R and full fusion.
    global alpr; global beta; global gamma; global sdelta; 
    global deltapr; global eta; global p; global qpr;
    global k; global f; global ep
    
    
    q=1-qpr; ppr=1-p;
    
    ep= lambd*((gamma**2)+(sdelta*gamma*ppr)+(lambd*ppr*sdelta))/(gamma*(gamma+lambd))
    print "Epsilon Value:\t%f" %(ep)
    
    theta= gamma +sdelta*(1+(lambd*p)/gamma)
    b_dot= theta - ep*(1+gamma/lambd)
    
    print "Theta:\t%f" %(theta)
    print "B_dot:\t%f" %(b_dot)
    
    r= ppr*sdelta*beta/(b_dot*(beta-k)-ppr*beta*sdelta)
    
    print "R:\t%f" %(r)
    
    num= k*((1+r)*deltapr*q*(b_dot+k-theta)+b_dot*(eta-k))
    #num= deltapr*q*(((r*b_dot)/(1+r))+k-theta) +eta-k
    #denom= (deltapr*q*gamma/k) -(1 - eta/k)*(b_dot/(1+r))
    denom=b_dot*(eta-k)-gamma*deltapr*q*(1+r)
    
    alpr=num/denom
    
    print "Alt Alpha Prime Value:\t%f" %(alpr)
    
def rootfinder3(): #Non-strict segregation of  K&R and full fusion.

    global alpr; global beta; global gamma; global sdelta; 
    global deltapr; global eta; global p; global qpr;
    global k; global f; global ep
    
    
    q=1-qpr; ppr=1-p;
    
    ep= lambd*((gamma**2)+(sdelta*gamma*ppr)+(lambd*ppr*sdelta))/(gamma*(gamma+lambd))
    print "Epsilon Value:\t%f" %(ep)
    
    theta= gamma +sdelta*(1+(lambd*p)/gamma)
    b_dot= theta - ep*(1+gamma/lambd)
    
    v= beta*((sdelta*ppr/b_dot)-1)+ ep-(p*sdelta*lambd)/gamma
    r= ppr*sdelta*beta/(b_dot*(beta-k)-ppr*beta*sdelta)
    
    tau= gamma*beta*( b_dot*(beta+v-k) - ppr*beta*sdelta )
    rho= (k-eta)*ppr*sdelta + gamma*r*deltapr*q
    phi=(b_dot/(ppr*sdelta))*(beta-k) - beta
    
    omg1=k*p*sdelta*b_dot*((lambd/gamma)+(beta/b_dot))
    omg2=(k*tau/rho)*(k-eta)*(1+((k-theta)/gamma))
    omg3=omg2 + ( k*tau*b_dot*(k-eta) )/(rho*gamma)
    omg4=(k*tau/rho)*( deltapr*q*(1+r) + (k-eta)*(b_dot/gamma) )
    
    #k_meu= (omg1-sdelta*beta)*(1)
    
    #num= k*((1+r)*deltapr*q*(b_dot+k-theta)+b_dot*(eta-k))
    #num= deltapr*q*(((r*b_dot)/(1+r))+k-theta) +eta-k
    #denom= (deltapr*q*gamma/k) -(1 - eta/k)*(b_dot/(1+r))
    #denom=b_dot*(eta-k)-gamma*deltapr*q*(1+r)
    
    alpr=(k/gamma)*( k -theta + b_dot*(1- (omg3+phi*(k-ep)*b_dot)/( phi*(omg1+b_dot*(k-ep)-sdelta*beta) +omg4 -v*beta*k*b_dot)) )
    
    print "Alt Alpha Prime Value:\t%f" %(alpr)
    
    #alpr=-0.0025
    
    a_dot=theta-k+(gamma*(alpr))/k
    
    k_meu_alpr= (omg1-sdelta*beta)*(1-a_dot/b_dot)+(ep-k)*a_dot
    
    alpha= (ppr*sdelta/tau)*( phi*k_meu_alpr + v*beta*k*b_dot*( (a_dot/b_dot) - 1))
    
    print "Alpha Value:\t%10.9f" %(alpha)
    
    alpha=(ppr*sdelta/tau)*(omg3 - omg4*(1-a_dot/b_dot))
    
    print "Alt Alpha Value:\t%10.9f" %(alpha)

init()
    
