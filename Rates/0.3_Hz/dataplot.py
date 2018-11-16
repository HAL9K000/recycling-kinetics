# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 17:40:21 2018

@author: Koustav
"""
import os
import numpy as np
from scipy.optimize import curve_fit as curvfit
#import math
import matplotlib.pyplot as plt


def f(x,a,b,c):
    return a+ b*np.exp(-c*x)

def main():
    gam=  0.06
    lam= 0.00405
    delt=1.0/120
    os.chdir("L_%f_g_%f_alpha_%f_d" %(lam, gam, delt))
    data=np.genfromtxt("Pre 0.3 Hz Alpha Rate Data.csv",comments="#",delimiter=',',skip_footer=2)
    xdata=data[:,0]
    ydata=data[:,4]
    plt.plot(xdata,ydata, marker='o', markerfacecolor='none', markeredgecolor='b', label="Points")
    popt, pcov = curvfit(f, xdata, ydata, bounds=([0.0,0,0],[0.75,0.1, 0.05]))
    print popt
    plt.plot(xdata, f(xdata, *popt), 'c-', label='fit: %5.4f + %5.4f*e^(-%5.4fx)' % tuple(popt))
    plt.xlabel("Time")
    plt.ylabel("dA/dt")
    plt.legend()
    plt.savefig("Pre Bestfit.png", dpi=200)
    plt.show()
    plt.close()


main()