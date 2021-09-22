# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 11:34:49 2020

@author: Boitshoko
""" 
import Fourier
import numpy as np
import matplotlib.pyplot as plt





n=500  #the number ofcells in the circuit
#constants
n=2*n+1
w=628318530         #frequency

Ri=50                #internal resistance
L=((270e-11)*1j*w)   #inductance
R=(30e-8)            #resistance
G=(1e7)              #ohmic losses
C=(1/(1j*w*110e-14))     #capacitance

R_L=L+R             #Series impedence
RL=100              #Terminal resistor
V_I=100
GC=1/(1/G+1/C)       #parallel impedance for C and G
Matrix=np.zeros((n,n), dtype=complex) #initialize the N*N MATRIX THAT will store impedandeces

def make(Matrix,w):
    #Make the matrix needed for calculations
    L=((270e-11)*1j*w)
    for i in range(1,n):
        
       if i%2==0 and i<n-1:
           Matrix[i][i]=L+R 
           Matrix[i][i+1]= 1/(1/G+1/C)
           Matrix[i+1][i+1]=-1
           Matrix[i+1][i+2]=-1
           Matrix[i][i-1]=-1/(1/G+1/C)
           Matrix[i+1][i]=1
    Matrix[0][0]=L+R+Ri
    Matrix[0][1]=1/(1/G+1/C)
    Matrix[1][2]=-1
    Matrix[n-1][n-2]=-1/(1/G+1/C)
    Matrix[n-1][n-1]=RL
    Matrix[1][0]=1
    Matrix[1][1]=-1
    return Matrix
Matrix=make(Matrix,w)

def Calculate(Matrix):
    
    '''This function solves for all currents in the wire 
       for a givin Voltage'''
    Voltage=np.zeros(n, dtype=complex)

    Voltage[0]=V_I
    I= np.linalg.solve(Matrix,Voltage)
    I_C=[]
    #I=np.absolute(I)
    for i in range(len(I)):
        if i%2==0:
            I_C.append(I[i])
    vol=[]
    vol2=[]
    for i in range(0,len(I_C)):
        vol.append((np.absolute(I_C[i])*R_L*(V_I*10/35)))
        vol2.append(I_C[i]*R_L*(V_I*10/35))
    return [np.absolute(vol),vol2]



def VoltagePm(Matrix):
    """This function calculats voltage along the the cable
       for different values of terminal resitor RL
       Quetion 2 b
    
    """

    #For open shorted cable set Terminal resistance RL to zero
    RL=0
    Matrix[n-1][n-1]=RL
    VS=Calculate(Matrix)[0]   #Stores voltages along the wire for shorted cable
    
    #For open ended cable set Terminal resistance RL to infinity

    RL=10e100
    Matrix[n-1][n-1]=RL
    VO=Calculate(Matrix)[0]    #Stores voltages along the wire for open cable

    #For open ended cable set Terminal resistance RL=100 Ohm
    RL=100
    Matrix[n-1][n-1]=RL
    VT=Calculate(Matrix)[0]    ##Stores voltages along the wire for 100 Ohm cable

    #Plot the results
    
    plt.plot(VO,label ='Open-Cable')
    plt.plot(VS,label ='Shorted-Cable')
    plt.plot(VT,label ="RL=100 Ohm")
    plt.legend(loc=2)
    plt.xlabel("x (cm)")
    plt.ylabel("Voltage (V)")
    plt.show()
    


#Plot voltage along the wire for the given matrix
VoltagePm(Matrix)

def RLmin(Matrix):
    '''The following code find the terminal resistace that minimizes the Amplitude
    of voltage along the transmision line
    '''
    #Stores minimum amplitudes of V for different values of termial resistance
    VMIN=[]  
    result=V_I
    index=-1    
    RTerm=[] #Stores different values of Terminal resitance
    k=0
    for i in range(40,60,1):
        """
        Here different values of RL are tryed
        """
        
        RL=i
        Matrix[n-1][n-1]=RL
        VT=Calculate(Matrix)[0]  
        a=np.amax(VT)-np.amin(VT)
        VMIN.append(a)
        RTerm.append(i)
        if(result>a):
            result=a
            index=k
        k+=1
    
    plt.plot(RTerm,VMIN)
    plt.xlabel("Terminal Resistance(Ohm)")
    plt.ylabel("Amplitude (V)")
    plt.show()
    
    return RTerm[index]
print(RLmin(Matrix))


def V_t(T,t,V_I):
    t0=3.3640559489972123e-09
    return V_I*np.exp((t0 - t)/T)*np.heaviside(t - t0, np.real(V_I))




N = 1000
t = np.linspace(0, 100e-9, N, endpoint=True)
V_t=V_t(5e-9,t,V_I)
Four_t=Fourier.dft(V_t)
freq = np.fft.fftfreq(t.shape[-1])
Vol=[]
for i in freq:
    Matrix=make(Matrix,2*np.pi*i*10e7)
    x=Calculate(Matrix)[1]  
    Vol.append(x[len(x)-1])
Vol=np.fft.ifft(Vol*1000)
    


plt.plot(t*1e9, V_t, label="Signal")
plt.xlabel("Time (ns)")
plt.ylabel('Voltage (V)')
plt.legend()
plt.show()

plt.plot(freq,Four_t)
plt.xlabel("Frequency (Hz)")
plt.show()


