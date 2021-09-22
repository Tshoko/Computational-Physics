# -*- coding: utf-8 -*-
"""
http://theflyingkeyboard.net/algorithms/python-trapezoidal-rule/
Created on Wed Apr 22 16:08:34 2020

@author: Boitshoko
"""
import scipy.integrate as integrate
#from scipy.integrate import quad
import json
from pylab import *
from scipy.special.orthogonal import p_roots
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np 
import scipy.integrate as spi
#Constants

a =-1; b = -a;

#f = lambda x : function(x)
#I = integrate.quad(f,a,b)
I=0.682689492137086
exact=I*b
# Calculates integral of f(x) from 0 to 1 for n
ns=[2,4,10,20,40,80,160]
#n for Gaussion quadrature method
ng=[2,3,4]

#Array to store Simpsons estimates values
simpsons=np.zeros(len(ns))
#Trapizoid estimates
trapizoid=np.zeros(len(ns))
#Gaussion quadrature estimates
Gaussion_Q=np.zeros(len(ng))
#Erros in Simspons estimates
err_S=np.zeros(len(ns))
#Erros in Troipizoid estimates
err_T=np.zeros(len(ns))
#Erros in Gaussion quadraturm estimates
err_G=np.zeros(len(ng))
# Richardson extrapolation for Simpsons and Traipizoid Method
rich_simps=np.zeros(len(ns))
rich_traps=np.zeros(len(ns))



def function(x):
    mu=0
    sigma=1
    
    y = scipy.stats.norm.pdf(x,mu,sigma)
    #y=(1/(sigma*np.sqrt(2*np.pi)))*(np.exp(-((x-mu)**2)/(2*sigma**2)))
    return y

def Richerdson_trap(A_n,A_2n):
    """Richerdson method for trapezodal rule"""
    return A_2n+(1/3)*(A_2n-A_n)

def Richerdson_Simps(A_n,A_2n):
    """Richerdson method for simpson rule"""

    return A_2n+(1/15)*(A_2n-A_n)  

def Trips(function,a,b,N):
    x=np.zeros(N)
    intervalBegin = a
    intervalEnd = b
    iterations = N
    step = (intervalEnd - intervalBegin) / iterations
    integral = 0.5 * (function(intervalBegin) + function(intervalEnd))
    for i in range(N):
        x[i]=intervalBegin + step * i
        
        integral += function(intervalBegin + step * i)
    integral *= step
    return integral



def Simpson(function, a, b, n):
    if a > b:
        print('Incorrect bounds')
        return None
    if n%2 != 0: # also an 'if' because both tests are NOT
        # mutually exclusive
        print('Invalid choice of n')
        return None
    else:
    	h = (b - a)/n # need to cast 'n' as float in order to avoid
        # integer division
    	sum1 = 0
    	for i in range(1, int(n/2 + 1)):
        	sum1 += function(a + (2*i - 1)*h)
    	sum1 *= 4
    	sum2 = 0
    	for i in range(1, int(n/2)): # range must be ints: range() integer 
            #end argument expected, got float.
        	sum2 += function(a + 2*i*h)
    	sum2 *= 2
    	approx = (b - a)/(3.0*n)*(function(a) + function(b) + sum1 + sum2)
    	return approx


i=0
for n in ns:
    
    s=Simpson(function, a, b, n)
    x = np.linspace(a,b,n+1)
    y = function(x)
    c = spi.simps(y,x)
    simpsons[i]=s
    err_S[i]=abs(c-exact)/exact
    t= np.trapz(y, x)
    trapizoid[i]=t
    err_T[i]=abs(t-exact)/exact
    i+=1  

def Richardson_ex():
    for j in range(0,len(ns)):
        if(j<len(ns)-1):
            rich_simps[j]=Richerdson_Simps(simpsons[j],simpsons[j+1])
            rich_traps[j]=Richerdson_trap(trapizoid[j],trapizoid[j+1])
        else:
        
            x = np.linspace(a,b,320+1)
            y = function(x)
            s = spi.simps(y,x)
            t= np.trapz(y, x)
            rich_simps[j]=Richerdson_Simps(simpsons[j],s)
            rich_traps[j]=Richerdson_trap(trapizoid[j],t)


def gauss(f,n,a,b):
    [x,w] = p_roots(n+1)
    
    G=0.5*(b-a)*sum(w*f(0.5*(b-a)*x+0.5*(b+a)))
    error=abs(G-exact)/exact
    return [G,error]
def Gauss(n):
    j=0
    for i in n:
        Gaussion_Q[j]=gauss(function,i,a,b)[0]
        err_G[j]=gauss(function,i,a,b)[1]
        j+=1
    return [Gaussion_Q,err_G]

def Json():
    Json={}
    nd=["2","4","10","20","40","80","160"]
    trap={}
    nd=ns
    sim={}
    RT={}
    RS={}
    Q={}
    for i in range(0,len(ns)):
        trap.update({ns[i]:trapizoid[i]})
        sim.update({nd[i]:simpsons[i]})
        #print(simpsons[i])
        RT.update({nd[i]:rich_traps[i]})
        RS.update({nd[i]:rich_simps[i]})
    for i in range(0,len(ng)):
        Q.update({str(ng[i]):Gaussion_Q[i]})
    a_1={"trapezoid":trap,"simpson":sim}
    b_1={"trapezoid":RT,"simpson":RS}
    
    Json.update({"1":{"a":a_1,"b":b_1,"c":Q},"2":{"a":"null"}})
    with open("tut3/MTSBOI022/tut3.json","w") as f:
        json.dump(Json,f)
    return Json




Richardson_ex()
Gauss(ng)
#create Json file
Json()


#plt.plot(ng,err_G,'*',label="Gaussian Quadrature")
plt.errorbar(ng,Gaussion_Q, yerr=err_G,fmt='o', label='Gaussian Quadrature and relative errors',)
plt.title("Relative error of numerical integration vs order")
plt.xlabel("Order (N) Of Integration")
plt.ylabel("Relative Error")
plt.legend()
plt.show()

plt.plot(ns,err_T,label="Trapezoidal Rule") 
plt.plot(ns,err_S,label="Simpson's Rule") 

plt.xlabel("Order (N)")
plt.ylabel("Relative Error")
plt.title("Relative Error Of Numerical Integration VS Order (N)")
plt.legend()
plt.show()    