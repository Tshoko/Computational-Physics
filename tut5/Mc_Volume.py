# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 20:59:15 2020

@author: boitshoko
"""

import matplotlib.pyplot as plt

import random
from scipy.special import gamma
import numpy as np
import time
import json

N=100
dim=3
Array=[]
def Vol(dim,N):
    """
      This Method uses mc  to calculate the
      volume for a given dimension and sample size.
    
    
    """
    hit=0
    
    for i in range(N):
        r=0
    
        for j in range(dim):
            # generate n random numbers
            #r = np.random.uniform(low=-1.0, high=1.0, size=(,))
            x=random.uniform(-1, 1)
    
            r += x*x
        if ((r)**(0.5) < 1):
            hit += 1
    return (2**(dim) * (hit /N))


def plot_mcint_stats(f,bins=10):
    """ This function takes a array of samples, interpreting the indices as runs and
      trials, respectively. It then plots some statistics """
    avg = np.average(f)
    #print('mean={',np.average(avg),'} ','std dev=',np.std(f),'}')
      
    #h,bins,p = plt.hist(f, bins=bins)
   
    # from now on we reuse the binning
    return [avg,np.std(f)/(1000)**(0.5)]

def Three_d_Volume(dim):
    """"
     This methpd give the estimation for volume
    
    """
    
    Array=[]
    for j in range(1000):
        v=Vol(dim,int(10))
        Array.append(v)
    return plot_mcint_stats(np.array(Array))
    '''
    plt.xlabel("estimate of Volume")
    plt.ylabel("#Samples")
    plt.show()'''

def Vol_Conver():
    """" 
    The Method Shows How the Volume estimations 
    Converges as the number of triels increase
    """
    l=np.linspace(10,500,100)
    for i in l:
        for j in range(10):
            v=Vol(dim,int(i))
            plt.plot(i,v,"*")
    plt.xlabel("#trials")
    plt.ylabel("Estimation of volume")
    plt.show()
def Time_estimate():
    times=[]
    means=[]
    std=[]
    true=[]
    for i in range (2,11):
        A=[]
        a=''
        print("Dimension",i)
        for k in range(100):
        
            start = time.time()
            
            a=Three_d_Volume(i)
    
            Volume=(np.pi**(i/2))/(gamma(i/2+1))
            end = time.time()
            
            elapsed = end - start
            A.append(elapsed)
        elapsed=np.average(np.array(A))
        times.append(elapsed)
        means.append(a[0])
        std.append(a[1])
        true.append(Volume)
        print("true Value is :",Volume)
        print("Estimated value is",a[0],"+/-",a[1])
        print("Avarage Computation time is:",elapsed)
    return [times,means,std,true]
      

Vol_Conver()
print(Three_d_Volume(3))
Array=Time_estimate()
x=[]
for i in range(2,11):
    x.append(i)
x=np.array(x)
print(x)
names = ["3d", "4d", "5d", "6d", "7d", "8d", "9d", "10d"]
output={}
k=3
v={}
u={}
for i in names:
    v[i]= Three_d_Volume(k)[0]
    u[i]= Three_d_Volume(k)[1]
    k+=1
out={
    "1": {
        "volume": v,
        "uncertainty": u
        
       }
   }
print(out)
def write_file():
    with open("./tut5/MTSBOI022/answers.json", 'w') as f:
        json.dump(out, f)
write_file()
#plot graphs
plt.plot(x,Array[0],"*-")
plt.xlabel("Dimension")
plt.ylabel("Avarage Computation time(seconds)")
plt.show()
#plot for volume
plt.plot(x,Array[1],"*-",label="Estimated Volume")
plt.plot(x,Array[3],"o-",label="Calculated Volume")
plt.legend()
plt.xlabel("Dimension")
plt.ylabel("Volume(Units)")
plt.show()
#plots for uncertainty
plt.plot(x,np.abs(np.array(Array[1])-np.array(Array[3])),"*-",label="fractional difference")

plt.plot(x,Array[2],"o-",label="Uncertainty of the mean")
plt.legend()
plt.xlabel("Dimension")
plt.ylabel("Error")


     