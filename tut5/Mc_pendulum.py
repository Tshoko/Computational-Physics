# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 23:28:05 2020

@author: boitshoko
Refference=https://nbviewer.jupyter.org/urls/www.numfys.net/media/notebooks/double_pendulum.ipynb
"""

from scipy.integrate import odeint

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

from numpy import cos, sin, arange, pi
import matplotlib.cm as cm
import numpy as np

#matplotlib inline
figsize = 6
dpi = 600


def RHS(t, z,Con ):
    L1,L2,m1,m2,g=Con
    theta1, theta2,w1, w2 = z
    cos12 = np.cos(theta1 - theta2)
    sin12 = np.sin(theta1 - theta2)
    sin1 = np.sin(theta1)
    sin2 = np.sin(theta2)
    xi = cos12**2*m2 - m1 - m2
    w1dot = ( L1*m2*cos12*sin12*w1**2 + L2*m2*sin12*w2**2
            - m2*g*cos12*sin2      + (m1 + m2)*g*sin1)/(L1*xi)
    w2dot = -( L2*m2*cos12*sin12*w2**2 + L1*(m1 + m2)*sin12*w1**2
            + (m1 + m2)*g*sin1*cos12  - (m1 + m2)*g*sin2 )/(L2*xi)
    return np.array([w1, w2, w1dot, w2dot])
def X(o1, o2,Con):
    l1,l2,m1,m2,g=Con
    x = l1*np.sin(o1) + l2*np.sin(o2)
    x1 = -l1*np.cos(o1) - l2*np.cos(o2)
    return [x,x1]
def to_cartesian(o1, o2,Con):
    L1,L2,m1,m2,g=Con
    """ Transforms theta and omega to cartesian coordinates
    and velocities x1, y1, x2, y2, vx1, vy1, vx2, vy2
    """
    theta1=o1
    theta2=o2
    x1 = L1 * np.sin(theta1)
    y1 = -L1 * np.cos(theta1)
    x2 = x1 + L2 * np.sin(theta2)
    y2 = y1 - L2 * np.cos(theta2)
    vx1 = L1*np.cos(theta1)*w1
    vy1 = L1*np.sin(theta1)*w1
    vx2 = vx1 + L2*np.cos(theta2)*w2
    vy2 = vy1 + L2*np.sin(theta2)*w2
    return [x1, y1, x2, y2, vx1, vy1, vx2, vy2]
def rk4(func,t,y0,Con):
    y = np.zeros((len(t),len(y0)))
    y[0] = y0
    x1=[]
    x2=[]
    x3=[]
    x4=[]
    d=[]
    for i in range(len(t)-1):
    
        dt = t[i+1] - t[i]
        
        k1 = dt * func( t[i]          , y[i]          ,Con)
        k2 = dt * func( t[i] + 0.5*dt , y[i] + 0.5*k1 ,Con)
        k3 = dt * func( t[i] + 0.5*dt , y[i] + 0.5*k2 ,Con)
        k4 = dt * func( t[i] + 1.0*dt , y[i] + 1.0*k1 ,Con)
        
        y[i+1] = y[i] + 1./6. * ( k1 + 2*k2 + 2*k3 + k4)
        d=to_cartesian(y[i][0],y[i][1],Con)
        x1.append(d[0])
        x2.append(d[1])
        x3.append(d[2])
        x4.append(d[3])
        

    return [y,np.array(x1),np.array(x2),np.array(x3),np.array(x4)]


def plot_phasespace(theta1, w1, theta2, w2):
    """ Creates a phase-space plot for the double pendulum
    for (theta, omega).
    
    theta1 : array-like, size(n,). The first angle in the double pendulum
    w1     : array-like, size(n,). Angular velocity of the first angle
    theta2 : array-like, size(n,). The second angle in the double pendulum
    w2     : array-like, size(n,). Angular velocity of the second angle
    """
    plt.title(r"Phase-space diagram, $\theta_{10}=%.1f$, $\theta_{20}=%.1f$ "%(theta1[0], theta2[0])
             + r"$\omega_{10}=%.1f$, $\omega_{20}=%.1f$"%(w1[0], w2[0]))
    plt.plot(theta1, w1, label=r"$i=1$")
    plt.plot(theta2, w2, label=r"$i=2$")
    plt.legend()
    plt.xlabel(r"$\theta_i$, [rad]")
    plt.ylabel(r"$\omega_i$, [rad/s]")
    xlim = [np.min(theta1), np.max(theta1), np.min(theta2), np.max(theta2)]
    plt.xlim(np.min(xlim), np.max(xlim))
    plt.show()
def plot_Time(theta1,theta2, t):
    """ Creates a phase-space plot for the double pendulum
    for (theta, omega).
    
    theta1 : array-like, size(n,). The first angle in the double pendulum
    w1     : array-like, size(n,). Angular velocity of the first angle
    theta2 : array-like, size(n,). The second angle in the double pendulum
    w2     : array-like, size(n,). Angular velocity of the second angle
    """
    plt.title(r"diagram for when, $\theta_{10}=%.1f$, $\theta_{20}=%.1f$ "%(theta1[0], theta2[0]))
    plt.plot(t,theta1, label=r"$i=1$")
    plt.plot(t, theta2, label=r"$i=2$")
    plt.legend()
    plt.xlabel(r"$time(seconds)$, [rad]")
    plt.ylabel(r"$\theta_i$, [rad]")
    
    plt.show()
def plot_position(x1, y1, x2, y2, theta1, theta2,w1,w2, t):
    """ Plots the motion of the double pendulum in the
    xy-plane, as well as the angles and the angular
    velocities as a function of time.
    
    x1     : array-like, size(n,). x-posision of mass 1
    y1     : array-like, size(n,). y-posision of mass 1
    x2     : array-like, size(n,). x-posision of mass 2
    y2     : array-like, size(n,). y-posision of mass 2
    theta1 : array-like, size(n,). The first angle in the double pendulum
    theta2 : array-like, size(n,). The second angle in the double pendulum
    t      : array-like, size(n,). Time
    """
    
    plt.figure(figsize=(2*figsize, figsize), dpi=dpi)
    # xy-plot
    L = 1.1*(l1 + l2)
    ax = plt.subplot(2, 2, (1, 3), autoscale_on=False, xlim=(-L, L), ylim=(-L, L))
    plt.title(r"posistion diagram, $\theta_{10}=%.1f$, $\theta_{20}=%.1f$ "%(theta1[0], theta2[0])
             + r"$\L1_{10}=%.1f$, $\omega_{20}=%.1f$"%(w1[0], w2[0]))
    ax.plot(x1, y1, label=r"Track $m_1$")
    ax.plot(x2, y2, label=r"Track $m_2$")
    ax.plot([0, x1[0], x2[0]], [0, y1[0], y2[0]], "-o", label="Initial position", c='k')
    plt.ylabel(r"$y/L$")
    plt.xlabel(r"$x/L$")
    ax.legend()
    plt.show()

l1=10
l2=1
m1=1
m2=10
g=9.8
theta1=np.pi
theta2=np.pi/3
w1=0
w2=0
Con=[l1,l2,m1,m2,9.8]
y0 = np.array([theta1,theta2, w1, w2])
t = np.arange(0, 50, 1e-2)
print("Solution calculated with",len(t),"steps")

y = rk4(RHS,t,y0,Con)

plt.plot(t,y[0][:,0])
plt.xlabel("Time(s)")
plt.ylabel(r"$\theta_1$, [rad]")
plt.show()
plt.plot(t,y[0][:,2])
plt.xlabel("Time(s)")
plt.ylabel(r"$\theta_2$, [rad]")
plt.show()
plot_phasespace(y[0][:,0], y[0][:,2], y[0][:,1], y[0][:,3])

plot_position(y[1], y[2], y[3], y[4],y[0][:,0],y[0][:,1],y[0][:,2], y[0][:,3], t)
plot_Time(y[0][:,0], y[0][:,1], t)
from scipy.fftpack import fft

xfrequencies = fft(y[2])
tfrequencies = np.linspace(0.0, 0.5*50, 50//2)
plt.figure()
plt.plot(tfrequencies, (2/50)*np.abs(xfrequencies[:50//2]))
plt.show()

