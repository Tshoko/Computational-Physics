# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:58:06 2020

@author: Boitshoko
"""

import numpy as np
import time
import No_numpy
import matplotlib.pyplot as plt
import sys
import argparse

from scipy.spatial import distance
L=[[5,10,20,50,100,500,1000],[5,10,20,50,100,500,700,800,900,1000,1300,1500],[5,10,20,50,100,500,700,800,900,1000,1300,1500,1700,1900,2000,2300]]
tnumpy=[]
tbrute=[]
def MKNumpy(n):
    M=np.zeros((n,n))
    x=np.arange(1, n+1, dtype=float)
    for i in range(1,n+1):
        for j in range(1,n+1):
            M[j-1][i-1]=((37.1*i+91.7*(j**2))%20.0)-10.0
    return [M,x]
def PLOT(tnumpy,tbrute,Lengths):
   
    plt.plot(Lengths,tnumpy,".")
    plt.plot(Lengths,tbrute,"*")
    plt.plot(Lengths,tnumpy,".",label='Numpy package')
    plt.plot(Lengths,tbrute,"*",label='Array implementation')
    plt.legend()
    plt.xlabel('Size n')
    plt.ylabel(' execution time(s)')
    plt.show()
    #plt.savefig('plot.png', dpi=300, bbox_inches='tight'))
def parabola(x, a, b, c):
    
    return a*x**2 + b*x + c
#for dotproduct
def x_tx(Lengths,tnumpy,tbrute):
    j=0
    for i in Lengths:
        #make matrix and vectorss
        M=MKNumpy(i)[0]
        x=MKNumpy(i)[1]
        M1=No_numpy.Make(i)[1]
        x1=No_numpy.Make(i)[0]
        start_time = time.time()
        x_t = x.T
        dot=x*x_t
        end_time = time.time()
        tnumpy.append(end_time-start_time)
        start_time1 = time.time()
        dot=No_numpy.Dot(x1,x1)
        end_time1 = time.time()
        tbrute.append(end_time1-start_time1)
        j+=1
    tnumpy=np.array(tnumpy, dtype='f')
    tbrute=np.array(tbrute, dtype='f')
    PLOT(tnumpy,tbrute,Lengths)
    
    
    
#print(M.dot(M))
def MX(Lengths,tnumpy,tbrute):
    j=0
    for i in Lengths:
        #make matrix and vectorss
        M=MKNumpy(i)[0]
        x=MKNumpy(i)[1]
        M1=No_numpy.Make(i)[1]
        x1=No_numpy.Make(i)[0]
        start_time = time.time()
        
        Mx=M*x
        end_time = time.time()
        tnumpy.append(end_time-start_time)
        start_time1 = time.time()
        Mx=No_numpy.Mul(M1,x1)
        end_time1 = time.time()
        tbrute.append(end_time1-start_time1)
        j+=1
    tnumpy=np.array(tnumpy, dtype='f')
    tbrute=np.array(tbrute, dtype='f')
    print(tnumpy)
    print(tbrute)
    PLOT(tnumpy,tbrute,Lengths)
    
def X_tMX(Lengths,tnumpy,tbrute):
    j=0
    for i in Lengths:
        #make matrix and vectorss
        M=MKNumpy(i)[0]
        x=MKNumpy(i)[1]
        M1=No_numpy.Make(i)[1]
        x1=No_numpy.Make(i)[0]
        start_time = time.time()
        
        xMx_mult= np.matmul(x,np.matmul(M, x.T))
        end_time = time.time()
        tnumpy.append(end_time-start_time)
        start_time1 = time.time()
        xMx_mult=No_numpy.XMX(M1,x1)
        end_time1 = time.time()
        tbrute.append(end_time1-start_time1)
        j+=1
    tnumpy=np.array(tnumpy, dtype='f')
    tbrute=np.array(tbrute, dtype='f')
    print(tnumpy)
    print(tbrute)
    PLOT(tnumpy,tbrute,Lengths)


def MM(Lengths,tnumpy,tbrute):
    j=0
    for i in Lengths:
        #make matrix and vectorss
        M=MKNumpy(i)[0]
        x=MKNumpy(i)[1]
        M1=No_numpy.Make(i)[1]
        x1=No_numpy.Make(i)[0]
        start_time = time.time()
        
        xMx_mult= np.matmul(M, M)
        end_time = time.time()
        tnumpy.append((end_time-start_time))
        start_time1 = time.time()
        xMx_mult=No_numpy.Multiply(M1)
        end_time1 = time.time()
        tbrute.append((end_time1-start_time1))
        j+=1
    tnumpy=np.array(tnumpy, dtype='f')
    tbrute=np.array(tbrute, dtype='f')
    print(tnumpy)
    print(tbrute)
    PLOT(tnumpy,tbrute,Lengths)
#depending on which graph or calculation you want comment out
print("This Program Plots Performance differences Of Matrix and vector Multiplications as n gets large. ")
print("The Perfomance difference is between numpy packages and matrix multiplication done using arrays and loops")
print()
print("Enter 1 for x_tx ")
print("Enter 2 for Mx ")
print("Enter 3 for x_tMx ")
print("Enter 4 for MM")
print()

A=eval(input("Enter oparation you want to Compare: "))
print("The are three diffent arrays to plot execution time vs n.bigger arrays take a long time to run")
print("Enter 0 for the shortest")
print("Enter 1 ")
print("Enter 2 for get more n data points ")
B=eval(input("Enter : "))
Lengths=L[B]
#A = int(sys.argv[1])
if (A==1):
    x_tx(Lengths,tnumpy,tbrute)
    print("excecuted")
elif (A==2):
    MX(Lengths,tnumpy,tbrute)
elif (A==3):
    X_tMX(Lengths,tnumpy,tbrute)
elif (A==4):
    X_tMX(Lengths,tnumpy,tbrute)
else:
    print("incorrect input")

#X_tMX(Lengths,tnumpy,tbrute)
#print(XMX(x_t,x,M))



