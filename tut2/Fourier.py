# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 17:02:03 2020

@author: Boitshoko
"""

import cmath
# Discrete fourier transform
def dft(x):
    t = []
    N = len(x)
    for k in range(N):
        a = 0
        for n in range(N):
            a += x[n]*cmath.exp(-2j*cmath.pi*k*n*(1/N))
        t.append(a)
    return t
# Inverse discrete fourier transform
def idft(t):
    x = []
    N = len(t)
    for n in range(N):
        a = 0
        for k in range(N):
            a += t[k]*cmath.exp(2j*cmath.pi*k*n*(1/N))
        a /= N
        x.append(a)
    return x