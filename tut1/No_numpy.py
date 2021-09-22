# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 00:21:47 2020

@author: Boitshoko
"""


def Make(n):
    MArr=[]
    x=[]
    
    for i in range(1,n+1):
        row=[]
        x.append(i+0.0)
        for j in range(1,n+1):
            Round=round(((((37.1*i+91.7*(j**2))%20.0)-10.0)),2)
            row.append(Round)
        MArr.append(row)
    return [x,MArr]
def Dot(x,x_t):
    s=0
    for i in range(len(x)):
       s+=x[i]*x_t[i]
       
    return round(s,2)
def Mul(M,x):
    row=[]
    
    for i in range(len(x)):
        k=0.0
        for j in range(len(x)):
                k+=M[j][i]*x[j]
        row.append(k)
       
    return row

def Multiply(M):
    n=len(M)
    res = [[0 for x in range(n)] for y in range(n)]  
  
    # explicit for loops 
    for i in range(len(M)): 
        for j in range(len(M[0])): 
            for k in range(len(M)): 
  
                # resulted matrix 
                res[j][i] += round(M[i][k] * M[k][j],2)
    return Round(res)

def Round(res):
    
    n=len(res)
    for i in range(0,n):
        for j in range(0,n):
            res[i][j]=round(res[i][j],2)
            
    return res

def XMX(M,x):
    return Dot(Mul(M,x),x)


'''
print(Mul(MArr,x))
print(Dot(Mul(MArr,x),x))
'''