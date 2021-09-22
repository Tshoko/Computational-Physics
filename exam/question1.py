# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 20:19:23 2020
@author: boitshoko
"""
import matplotlib.pyplot as plt


from numpy import cos, sin, arange, pi
import matplotlib.cm as cm
import numpy as np
import json



def RHS(y,c):
    #general velocity
    g0=y[3]
    g1=y[4]
    g2=y[5]
    
    C=1
    K=223
    e=1
    ptot=(g0**2+g1**2+g2**2)
    mass=935
    c=(K*e*C**2)/((mass**2)*(C**4)+ptot*(C**2))**(0.5)
    
    g3=-c*y[5]*y[0]
    g4=c*y[5]*y[1]
    g5=c*(y[0]*y[3]-y[4]*y[1])
    

    return np.array([g0,g1,g2,g3,g4,g5])
def rk4(func,t,y0,Con):
    #4th order runger Kutter
    y = np.zeros((len(t),len(y0)))
    y[0] = y0
    x1=[y0[0]]
    x2=[y0[1]]
    x3=[y0[2]]
    x4=[y0[3]]
    x5=[y0[4]]
    x6=[y0[5]]
    d=[]
    for i in range(len(t)-1):
    
        dt = t[i+1] - t[i]
        
        k1 = dt * RHS( y[i]          ,Con)
        k2 = dt * RHS( y[i] + 0.5*k1 ,Con)
        k3 = dt * RHS( y[i] + 0.5*k2 ,Con)
        k4 = dt * RHS( y[i] + 1.0*k1 ,Con)
        
        y[i+1] = y[i] + 1./6. * ( k1 + 2*k2 + 2*k3 + k4)
        d=y[i]
        x1.append(d[0])
        x2.append(d[1])
        x3.append(d[2])
        x4.append(d[3])
        x5.append(d[4])
        x6.append(d[5])
        if(d[2]>3.1):
            break
        
    return [y,np.array(x1),np.array(x2),np.array(x3),np.array(x4),np.array(x5),np.array(x6)]
def Ploty(y0,cord):
    #plot trajectories
    mom='Momentum_'+cord
    pos=cord+'('+'m'+')'
    t = np.arange(0, 10, 1e-4 )
    print("Solution calculated with",len(t),"steps")
    c=1
    y = rk4(RHS,t,y0,c)
    print(y[1],y[4],y[6])
    plt.plot(y[3],y[2],label="Position")

    plt.title('Particles Trajectory') 
    plt.xlabel('distance Z(m)') 
    plt.ylabel(pos) 
    plt.show()
    plt.plot(y[3],y[5],label="Position")
    #Momentum
    plt.title('Particles Momentum') 
    plt.xlabel('distance Z(m)') 
    plt.ylabel(mom) 
    plt.show()
def Plot(y0,cord):
    #plot trajectories
    mom='Momentum_'+cord
    pos=cord+'('+'m'+')'
    t = np.arange(0, 10, 1e-4 )
    print("Solution calculated with",len(t),"steps")
    c=1
    y = rk4(RHS,t,y0,c)
    print(y[1],y[4],y[6])
    plt.plot(y[3],y[1],label="Position")

    plt.title('Particles Trajectory') 
    plt.xlabel('distance Z(m)') 
    plt.ylabel(pos) 
    plt.show()
    plt.plot(y[3],y[4],label="Position")
    #Momentum
    plt.title('Particles Momentum') 
    plt.xlabel('distance Z(m)') 
    plt.ylabel(mom) 
    plt.show()

    
    
""""Question B"""
cord='x'
Pz=1
y0 = np.array([0.01,0,0,0,0,Pz])
Plot(y0,cord)
Pz=1
y0 = np.array([-0.01,0,0,0,0,Pz])
Plot(y0,cord)
""""Question c a"""
cord='x'
Pz=1
y0 = np.array([0.1,0,0,0,0,Pz])
Plot(y0,cord)
Pz=1
y0 = np.array([-0.1,0,0,0,0,Pz])
Plot(y0,cord)
""""Question c b"""
print("Question c b")
cord='y'
Pz=1
y0 = np.array([0,-0.1,0,0,0,Pz])
Ploty(y0,cord)
Pz=1
y0 = np.array([0,0.1,0,0,0,Pz])
Ploty(y0,cord)
""""Question c b"""
print("Question c b")
cord='y'
Pz=1
y0 = np.array([0,-0.1,0,0,0,Pz])
Ploty(y0,cord)
Pz=1
y0 = np.array([0,0.1,0,0,0,Pz])
Ploty(y0,cord)

y0 = np.array([0,0,1,-0.71,0,0.71])
Plot(y0,'x')

JSON={
  "1": {
    "b": {
      "x": 0,
      "y": 0,
      "z": 1.55,
      "px": 0,
      "py": 0,
      "pz": 7000
    },
    "e": {
      "Bx": 0,
      "By": 0,
      "Bz": 0
    }
  }
}

def write_file():
    """ this writes to the the json file"""
    with open("./exam/MTSBOI022/answers.json", 'w') as f:
        json.dump(JSON, f)
write_file()