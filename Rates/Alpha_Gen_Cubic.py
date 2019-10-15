# -*- coding: utf-8 -*-
"""
Created on Mon Sep 09 21:00:11 2019

@author: Koustav
"""

import cmath
import os
import math
import numpy as np

alpr=alph=beta=gamma=delt=k=lambd=ep=0.0
p=0.5
ppr=1.0-p
print "Enter Frequency Number: \t",
f=float(raw_input())
if( os.path.isdir("%2.1f_Hz" %(f))==False):
    os.mkdir("%2.1f_Hz" %(f))
os.chdir("%2.1f_Hz" %(f))

def main():
    
    global alpr; global alph; global beta; global ep; global gamma; global delt; global k; global f
    
    global p; global ppr
    
    #alpr=float(raw_input("Enter Alpha Prime Value: \t"))
    alpr= 0.0025
    #beta=float(raw_input("Enter Beta Value: \t"))
    beta=1.0/5.0
    lambd=0.0621
    gamma=0.2*f
    delt=1.0/120.0
    ep=(lambd/(gamma*(gamma+ lambd)))*((gamma)**2 + (ppr*delt)*(gamma+lambd))
    #delta= float(raw_input("Enter Undocking Rate Value: \t"))
    
    k= float(raw_input("Enter K Value: \t"))
    
    m= 1 + gamma/k
    theta= delt*(1+ppr*lambd/gamma)
    b_pr= gamma + theta -ep*(1 + gamma/k)
    omeg= 1 + lambd/gamma -ep- p*(lambd/gamma + delt/b_pr)
    
    n= alpr*(gamma+theta+m*(1-alpr))+ delt*p+ b_pr*lambd*ppr/gamma -(gamma+theta)
    
    print "omega:\t%5.3f\tn:\t%5.3f\tb_pr:\t%5.3f\tm:\t%5.3f\t" %(omeg,n,b_pr,m)
    
    #Coefficients of the cubic polynomial
    
    a=1
    b= (2+m)*alpr -(m+gamma+theta+ep)
    c= (gamma+theta+m - (1+m)*alpr)*(ep-alpr) - n
    d= n*(ep-alpr)- lambd*ppr*delt*b_pr*omeg/gamma
    
    print "b:\t%5.3f\tc:\t%5.3f\td:\t%5.3f\t" %(b,c,d)
    
    '''Applying Tschinachius' Transformation here'''
    
    p= (3*a*c-b**2)/(3*a*a)
    q= (2*b**3-9*a*b*c + 27*(a*a)*d)/(27*a**3)
    
    chec1= 4*p**3 +27*q**2
    chec2= 18*a*b*c*d - 4*(b**3)*d+ (b**2)*(c**2)- 4*a*(c**3) -27*(a*d)**2
    
    print("p q check:\t%4f" %(chec1))
    print("Discriminant check:\t%f" %(chec2))
    
    i= math.pow((-((q/2)**2 +(p/3)**3)), 0.5)
    #Cube root is imaginary
    
    z= complex(-q/2,i)
    z_con= z.conjugate()
    
    za=pow(z, 1/3)
    
    print 'Predicted Cube Root Of Z:\t'+str(za)
    
    print "Z:\t"+str(z)+"\tZ Conj:\t"+str(z_con)
    
    r,phi = cmath.polar(z)
    r_con, phi_con = cmath.polar(z_con)
    #Converting To Polar Co-ordinates
    
    r= math.pow(r, 1/3)
    r_con=math.pow(r_con, 1/3)
    
    phi1=[]; phi_con1=[]
    
    for x in range(-2,3,2):
        print"x:\t%d" %(x)
        phi1.append((phi/3+x*cmath.pi/3))
        phi_con1.append((phi_con/3+x*cmath.pi/3))
    #Generating the phase values of all the different cube roots.
    
    flag=0; flag_con=0  #Stores maximum real value
    
    print "Possible Cubic Roots Of Z"
    
    for t in range(0,len(phi1)):
        
        z1=cmath.rect(r, phi1[t])
        
        print z1
        
        z1_con=cmath.rect(r_con, phi_con1[t])
        
        if z1.real>= flag:
            flag=z1.real
            z=z1
        
        
        if z1_con.real>= flag_con:
            flag_con=z1_con.real
            z_con=z1_con
            
        
    print "Maximal Principal Cube Root Of"
    print "Z:\t"+str(z)
    print "Z Conjugate:\t"+str(z_con)
            
    print "Possible Roots of Depressed Cubic"
    
    omg1=1
    omg2=complex(-1/2, (math.pow(3,1/2)/2))
    omg3=complex(-1/2, -(math.pow(3,1/2)/2))
    
    #Cube Roots of Unity
    
    t=[]
    #Stores the roots of depressed cubic
    
    r=z_con+z
    t.append(r.real-b/(3*a))
    print r.real
    
    r=omg2*z +omg3*z_con
    t.append(r.real-b/(3*a))
    print r.real
    
    r=omg3*z +omg2*z_con
    t.append(r.real-b/(3*a))
    print r.real
    
    print "Possible Roots of Actual Cubic"
    
    print t
    
    print "Checking by Other Method"
    
    roots= np.roots([a,b,c,d])
    
    print roots
    
    
    
    
    
    
        
    
    
main()