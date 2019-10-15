# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 16:08:59 2019

@author: Koustav
"""

import os
import numpy as np


class Multi_Alpha_Gen_Cubic:
    
    
    def __init__(self):
    
        self.p=0.5 
        self.ppr=0.5
    
    
    
        self.f=int(raw_input("Enter the input frquency ( must be either 1, 3, 10 or 30):\t"))
        #self.alpr=float(raw_input("Enter Alpha Prime Value: \t"))
        self.alpr= 0.0025
        #self.beta=float(raw_input("Enter self.beta Value: \t"))
        self.beta=1.0/5.0
        self.lambd={10: 0.0310, 30: 0.0621, 3: 0.0132, 1: 0.0056}
        self.gamma=0.2*self.f
        self.delta=1.0/120.0
        self.ep=(self.lambd[self.f]/(self.gamma*(self.gamma+ self.lambd[self.f])))*((self.gamma)**2 + (self.ppr*self.delta)*(self.gamma+self.lambd[self.f]))
        #self.delta= float(raw_input("Enter Undocking Rate Value: \t"))
    
        self.k= {10: 0.01489, 30: 0.0538, 3: 0.00571, 1: 0.00456}
        
        self.roots=[] #Stores data on roots.
    
    
    def rootfinder(self):
        
        j=0  #Used to update self.p
        
        i=0 #Used as a simple counter to determine the number of times the loop below has executrd.
        
        #Change to relevant directory
        if( os.path.isdir("%2.1f_Hz" %(self.f))==False):
            os.mkdir("%2.1f_Hz" %(self.f))
        os.chdir("%2.1f_Hz" %(self.f))
        
        if(os.path.isdir("Varying p")==False):
            os.mkdir("Varying p")
        os.chdir("Varying p")
        
        
        self.data_exp=np.zeros([11,4])
        #Array used to export data generated on roots
            
        
        while j<=1.0:
            
            self.p=j; self.ppr= 1-j
            self.ep=(self.lambd[self.f]/(self.gamma*(self.gamma+ self.lambd[self.f])))*((self.gamma)**2 + (self.ppr*self.delta)*(self.gamma+self.lambd[self.f]))
    
            m= 1 + self.gamma/self.k[self.f]
            theta= self.delta*(1+self.ppr*self.lambd[self.f]/self.gamma)
            b_pr= self.gamma + theta -self.ep*(1 + self.gamma/self.k[self.f])
            omeg= 1 + self.lambd[self.f]/self.gamma -self.ep- self.p*(self.lambd[self.f]/self.gamma + self.delta/b_pr)
    
            n= self.alpr*(self.gamma+theta+m*(1-self.alpr))+ self.delta*self.p+ b_pr*self.lambd[self.f]*self.ppr/self.gamma -(self.gamma+theta)
    
            print "omega:\t%5.3f\tn:\t%5.3f\tb_pr:\t%5.3f\tm:\t%5.3f\t" %(omeg,n,b_pr,m)
    
            #Coefficients of the cubic polynomial
    
            a=1
            b= (2+m)*self.alpr -(m+self.gamma+theta+self.ep)
            c= (self.gamma+theta+m - (1+m)*self.alpr)*(self.ep-self.alpr) - n
            d= n*(self.ep-self.alpr)- self.lambd[self.f]*self.ppr*self.delta*b_pr*omeg/self.gamma
    
            print "b:\t%5.3f\tc:\t%5.3f\td:\t%5.3f\t" %(b,c,d)
    
            '''Applying Tschinachius' Transformation here'''
    
            p= (3*a*c-b**2)/(3*a*a)
            q= (2*b**3-9*a*b*c + 27*(a*a)*d)/(27*a**3)
    
            chec1= 4*p**3 +27*q**2
            chec2= 18*a*b*c*d - 4*(b**3)*d+ (b**2)*(c**2)- 4*a*(c**3) -27*(a*d)**2
    
            print("p q check:\t%4f" %(chec1))
            print("Discriminant check:\t%f" %(chec2))
    
            
            print "P-value:\t%2.1f" %(self.p)
            print "Checking by Other Method"
    
            roots= np.roots([a,b,c,d])
    
            print roots
            
            self.roots.append(roots)
            
            self.data_exp[i,0]= self.p
            #First col of self.data_exp stores p values. 
            
            i+=1
            j+=0.1
        
        self.roots=np.array(self.roots)
        #Converting to numpy format for easy indexing
        
        for j in range(1,4):
            self.data_exp[:,j]=self.roots[:,j-1]
        
        #Here the three cubic roots obtained in each iteration get assigned to the 1st, 2nd and 3rd column of data_exp 1,2,3
        
        np.savetxt('Cubic Root Information.csv', self.data_exp, delimiter=",", header="0: p value, 1: 1st Root, 2: 2nd Root, 3: 3rd Root",comments="#")
            
        
        
obj=Multi_Alpha_Gen_Cubic()

if __name__ == '__main__':
    obj.rootfinder()
      
        
        